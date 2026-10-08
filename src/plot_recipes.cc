#include "plot_recipes.h"
#include "histograms.h"
#include <memory>
#include <vector>
#include <string>

namespace {

// Categories that are also plotted summed together, on top of the 13
// individual ones. These are the pairs whose distributions are close enough
// in shape -- and thin enough in statistics on their own -- that the cuts in
// pt2_cuts.h tune them as a single group. Plotting them the same way means
// the plots show what the cuts actually act on.
//
// To combine a different set, edit this list: each inner list is one plot.
const std::vector<std::vector<int>> kMergedGroups = {
    {0, 2},   // L1F_to_L2F + L1T_to_L2F
    {1, 3},   // L1F_to_L2T + L1T_to_L2T
    {5, 7},   // L2F_to_L3F + L2T_to_L3F
    {6, 8},   // L2F_to_L3T + L2T_to_L3T
};

// A const HistogramManager hands out its [13][2] grids as this.
using ConstHistGrid = TH1D *const (*)[kNCharge];

// Summed histograms made for the merged plots. PlotRecipe only borrows its
// pointers, so the sums have to outlive getPt2Recipes(); they live here until
// the process exits.
std::vector<std::unique_ptr<TH1D>> &mergedStore()
{
    static std::vector<std::unique_ptr<TH1D>> store;
    return store;
}

// The histogram to plot for one group: the original when the group names a
// single (category, charge), otherwise a fresh sum over every member.
TH1D *grouped(ConstHistGrid grid, const std::vector<int> &cats, const std::vector<int> &charges)
{
    if (cats.size() == 1 && charges.size() == 1) return grid[cats[0]][charges[0]];

    auto *sum = static_cast<TH1D *>(grid[cats[0]][charges[0]]->Clone());
    sum->SetDirectory(nullptr);
    sum->Reset();
    for (int i : cats)
        for (int q : charges) sum->Add(grid[i][q]);
    mergedStore().emplace_back(sum);
    return sum;
}

} // namespace

std::vector<PlotRecipe> getPt2Recipes(const HistogramManager& hists) {
    std::vector<PlotRecipe> recipes;

    // One entry per plot: which categories and charges to add together, and
    // how to name the result.
    struct PlotGroup {
        std::vector<int> cats;
        std::vector<int> charges;
        std::string sfx;
        std::string ttl;
    };
    std::vector<PlotGroup> groups;

    auto addGroup = [&](const std::vector<int> &cats, int c) {
        std::string names = hists.catNames[cats[0]];
        std::string titles = hists.catTitles[cats[0]];
        for (size_t k = 1; k < cats.size(); ++k) {
            names  += "_and_" + hists.catNames[cats[k]];
            titles += " + "   + hists.catTitles[cats[k]];
        }
        groups.push_back({cats, {c},
                          "_" + names + "_" + hists.chargeNames[c],
                          " (" + titles + ", " + hists.chargeTitles[c] + ")"});
    };

    // The 13 categories on their own, then each merged group.
    for (int i = 0; i < 13; ++i)
        for (int c = 0; c < 2; ++c) addGroup({i}, c);
    for (const auto &m : kMergedGroups)
        for (int c = 0; c < 2; ++c) addGroup(m, c);

    // Everything pooled: all 13 categories and both charges summed into one
    // plot per variable. The per-category pages show where a variable
    // separates real from fake; this one shows what the sample as a whole
    // looks like, which is what gets handed downstream.
    std::vector<int> allCats(13);
    for (int i = 0; i < 13; ++i) allCats[i] = i;
    groups.push_back({allCats, {0, 1}, "_AllCategories",
                      " (all categories, both charges)"});

    for (const auto &g : groups) {
        const std::string &sfx = g.sfx;
        const std::string &ttl = g.ttl;

        // =====================================================================
        // ALL pT2s - KINEMATICS & LST VARIABLES
        // =====================================================================

        recipes.push_back({
            .title = "pT2 Delta p_{T}" + ttl,
            .xAxis = "#Delta p_{T} [GeV]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPT" + sfx,
            .hists = {grouped(hists.real_pt2_deltaPT, g.cats, g.charges), grouped(hists.fake_pt2_deltaPT, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPHI" + sfx,
            .hists = {grouped(hists.real_pt2_deltaPHI, g.cats, g.charges), grouped(hists.fake_pt2_deltaPHI, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta #eta" + ttl,
            .xAxis = "#eta",
            .yAxis = "Entries",
            .filename = "pt2_all_ETA" + sfx,
            .hists = {grouped(hists.real_pt2_deltaETA, g.cats, g.charges), grouped(hists.fake_pt2_deltaETA, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // --- pLS Absolute Eta ---
        recipes.push_back({
            .title = "pLS Absolute #eta" + ttl,
            .xAxis = "pLS #eta",
            .yAxis = "Entries",
            .filename = "pt2_all_pls_ETA" + sfx,
            .hists = {grouped(hists.real_pt2_pls_ETA, g.cats, g.charges), grouped(hists.fake_pt2_pls_ETA, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pLS Absolute #eta" + ttl,
            .xAxis = "pLS #eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_pls_ETA" + sfx,
            .hists = {grouped(hists.real_unused_pt2_pls_ETA, g.cats, g.charges), grouped(hists.fake_unused_pt2_pls_ETA, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // --- LS Absolute Eta ---
        recipes.push_back({
            .title = "LS Absolute #eta" + ttl,
            .xAxis = "LS #eta",
            .yAxis = "Entries",
            .filename = "pt2_all_ls_ETA" + sfx,
            .hists = {grouped(hists.real_pt2_ls_ETA, g.cats, g.charges), grouped(hists.fake_pt2_ls_ETA, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LS Absolute #eta" + ttl,
            .xAxis = "LS #eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_ls_ETA" + sfx,
            .hists = {grouped(hists.real_unused_pt2_ls_ETA, g.cats, g.charges), grouped(hists.fake_unused_pt2_ls_ETA, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta R" + ttl,
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaR" + sfx,
            .hists = {grouped(hists.real_pt2_deltaR, g.cats, g.charges), grouped(hists.fake_pt2_deltaR, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdBeta" + sfx,
            .hists = {grouped(hists.real_pt2_LSTdBeta, g.cats, g.charges), grouped(hists.fake_pt2_LSTdBeta, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTKinZRes" + sfx,
            .hists = {grouped(hists.real_pt2_LSTKinZRes, g.cats, g.charges), grouped(hists.fake_pt2_LSTKinZRes, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTOrgZRes" + sfx,
            .hists = {grouped(hists.real_pt2_LSTOrgZRes, g.cats, g.charges), grouped(hists.fake_pt2_LSTOrgZRes, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST #Delta#phi" + ttl,
            .xAxis = "LST #Delta#phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdPhi" + sfx,
            .hists = {grouped(hists.real_pt2_LSTdPhi, g.cats, g.charges), grouped(hists.fake_pt2_LSTdPhi, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
         });


        // =====================================================================
        // ALL pT2s - MD COMPONENTS
        // =====================================================================

        recipes.push_back({
            .title = "MD0 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dxy" + sfx,
            .hists = {grouped(hists.real_pt2_MD0_dXY, g.cats, g.charges), grouped(hists.fake_pt2_MD0_dXY, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dz" + sfx,
            .hists = {grouped(hists.real_pt2_MD0_dZ, g.cats, g.charges), grouped(hists.fake_pt2_MD0_dZ, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dxy" + sfx,
            .hists = {grouped(hists.real_pt2_MD1_dXY, g.cats, g.charges), grouped(hists.fake_pt2_MD1_dXY, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dz" + sfx,
            .hists = {grouped(hists.real_pt2_MD1_dZ, g.cats, g.charges), grouped(hists.fake_pt2_MD1_dZ, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_rz_simple" + sfx,
            .hists = {grouped(hists.real_pt2_MD0_rz_simple, g.cats, g.charges), grouped(hists.fake_pt2_MD0_rz_simple, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_rz_simple" + sfx,
            .hists = {grouped(hists.real_pt2_MD1_rz_simple, g.cats, g.charges), grouped(hists.fake_pt2_MD1_rz_simple, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // =====================================================================
        // UNUSED pT2s - KINEMATICS & LST VARIABLES
        // =====================================================================

        recipes.push_back({
            .title = "Unused pT2 Delta p_{T}" + ttl,
            .xAxis = "#Delta p_{T} [GeV]",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPT" + sfx,
            .hists = {grouped(hists.real_unused_pt2_deltaPT, g.cats, g.charges), grouped(hists.fake_unused_pt2_deltaPT, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPHI" + sfx,
            .hists = {grouped(hists.real_unused_pt2_deltaPHI, g.cats, g.charges), grouped(hists.fake_unused_pt2_deltaPHI, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Absolute #eta" + ttl,
            .xAxis = "#eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaETA" + sfx,
            .hists = {grouped(hists.real_unused_pt2_deltaETA, g.cats, g.charges), grouped(hists.fake_unused_pt2_deltaETA, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Delta R" + ttl,
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaR" + sfx,
            .hists = {grouped(hists.real_unused_pt2_deltaR, g.cats, g.charges), grouped(hists.fake_unused_pt2_deltaR, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
           .printYields = true
        });


        recipes.push_back({
            .title = "Unused LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdBeta" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTdBeta, g.cats, g.charges), grouped(hists.fake_unused_pt2_LSTdBeta, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTKinZRes" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTKinZRes, g.cats, g.charges), grouped(hists.fake_unused_pt2_LSTKinZRes, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTOrgZRes" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTOrgZRes, g.cats, g.charges), grouped(hists.fake_unused_pt2_LSTOrgZRes, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST #Delta#phi" + ttl,
            .xAxis = "LST #Delta#phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdPhi" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTdPhi, g.cats, g.charges), grouped(hists.fake_unused_pt2_LSTdPhi, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // =====================================================================
        // UNUSED pT2s - MD COMPONENTS
        // =====================================================================

        recipes.push_back({
            .title = "Unused MD0 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dxy" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD0_dXY, g.cats, g.charges), grouped(hists.fake_unused_pt2_MD0_dXY, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dz" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD0_dZ, g.cats, g.charges), grouped(hists.fake_unused_pt2_MD0_dZ, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dxy" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD1_dXY, g.cats, g.charges), grouped(hists.fake_unused_pt2_MD1_dXY, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dz" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD1_dZ, g.cats, g.charges), grouped(hists.fake_unused_pt2_MD1_dZ, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_rz_simple" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD0_rz_simple, g.cats, g.charges), grouped(hists.fake_unused_pt2_MD0_rz_simple, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_rz_simple" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD1_rz_simple, g.cats, g.charges), grouped(hists.fake_unused_pt2_MD1_rz_simple, g.cats, g.charges)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });
    }

    return recipes;
}
