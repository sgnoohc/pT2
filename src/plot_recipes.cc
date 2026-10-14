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

// The histogram to plot for one group: the original when the group is a
// single category, otherwise a fresh sum over its members.
TH1D *grouped(ConstHistGrid grid, const std::vector<int> &cats, int c)
{
    if (cats.size() == 1) return grid[cats[0]][c];

    auto *sum = static_cast<TH1D *>(grid[cats[0]][c]->Clone());
    sum->SetDirectory(nullptr);
    for (size_t k = 1; k < cats.size(); ++k) sum->Add(grid[cats[k]][c]);
    mergedStore().emplace_back(sum);
    return sum;
}

} // namespace

std::vector<PlotRecipe> getPt2Recipes(const HistogramManager& hists) {
    std::vector<PlotRecipe> recipes;

    // The 13 categories on their own, then each merged group.
    std::vector<std::vector<int>> groups;
    for (int i = 0; i < 13; ++i) groups.push_back({i});
    for (const auto &g : kMergedGroups) groups.push_back(g);

    for (const auto &grp : groups) {
        for (int c = 0; c < 2; ++c) {
        std::string cats = hists.catNames[grp[0]];
        std::string titles = hists.catTitles[grp[0]];
        for (size_t k = 1; k < grp.size(); ++k) {
            cats += "_and_" + hists.catNames[grp[k]];
            titles += " + " + hists.catTitles[grp[k]];
        }
        std::string sfx = "_" + cats + "_" + hists.chargeNames[c];
        std::string ttl = " (" + titles +  ", " + hists.chargeTitles[c] + ")";

        // =====================================================================
        // ALL pT2s - KINEMATICS & LST VARIABLES
        // =====================================================================

        recipes.push_back({
            .title = "pT2 Delta p_{T}" + ttl,
            .xAxis = "#Delta p_{T} [GeV]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPT" + sfx,
            .hists = {grouped(hists.real_pt2_deltaPT, grp, c), grouped(hists.fake_pt2_deltaPT, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPHI" + sfx,
            .hists = {grouped(hists.real_pt2_deltaPHI, grp, c), grouped(hists.fake_pt2_deltaPHI, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta #eta" + ttl,
            .xAxis = "#eta",
            .yAxis = "Entries",
            .filename = "pt2_all_ETA" + sfx,
            .hists = {grouped(hists.real_pt2_deltaETA, grp, c), grouped(hists.fake_pt2_deltaETA, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // --- pLS Absolute Eta ---
        recipes.push_back({
            .title = "pLS Absolute #eta" + ttl,
            .xAxis = "pLS #eta",
            .yAxis = "Entries",
            .filename = "pt2_all_pls_ETA" + sfx,
            .hists = {grouped(hists.real_pt2_pls_ETA, grp, c), grouped(hists.fake_pt2_pls_ETA, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pLS Absolute #eta" + ttl,
            .xAxis = "pLS #eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_pls_ETA" + sfx,
            .hists = {grouped(hists.real_unused_pt2_pls_ETA, grp, c), grouped(hists.fake_unused_pt2_pls_ETA, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // --- LS Absolute Eta ---
        recipes.push_back({
            .title = "LS Absolute #eta" + ttl,
            .xAxis = "LS #eta",
            .yAxis = "Entries",
            .filename = "pt2_all_ls_ETA" + sfx,
            .hists = {grouped(hists.real_pt2_ls_ETA, grp, c), grouped(hists.fake_pt2_ls_ETA, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LS Absolute #eta" + ttl,
            .xAxis = "LS #eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_ls_ETA" + sfx,
            .hists = {grouped(hists.real_unused_pt2_ls_ETA, grp, c), grouped(hists.fake_unused_pt2_ls_ETA, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta R" + ttl,
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaR" + sfx,
            .hists = {grouped(hists.real_pt2_deltaR, grp, c), grouped(hists.fake_pt2_deltaR, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdBeta" + sfx,
            .hists = {grouped(hists.real_pt2_LSTdBeta, grp, c), grouped(hists.fake_pt2_LSTdBeta, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTKinZRes" + sfx,
            .hists = {grouped(hists.real_pt2_LSTKinZRes, grp, c), grouped(hists.fake_pt2_LSTKinZRes, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTOrgZRes" + sfx,
            .hists = {grouped(hists.real_pt2_LSTOrgZRes, grp, c), grouped(hists.fake_pt2_LSTOrgZRes, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST #Delta#phi" + ttl,
            .xAxis = "LST #Delta#phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdPhi" + sfx,
            .hists = {grouped(hists.real_pt2_LSTdPhi, grp, c), grouped(hists.fake_pt2_LSTdPhi, grp, c)},
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
            .hists = {grouped(hists.real_pt2_MD0_dXY, grp, c), grouped(hists.fake_pt2_MD0_dXY, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dz" + sfx,
            .hists = {grouped(hists.real_pt2_MD0_dZ, grp, c), grouped(hists.fake_pt2_MD0_dZ, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dxy" + sfx,
            .hists = {grouped(hists.real_pt2_MD1_dXY, grp, c), grouped(hists.fake_pt2_MD1_dXY, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dz" + sfx,
            .hists = {grouped(hists.real_pt2_MD1_dZ, grp, c), grouped(hists.fake_pt2_MD1_dZ, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_rz_simple" + sfx,
            .hists = {grouped(hists.real_pt2_MD0_rz_simple, grp, c), grouped(hists.fake_pt2_MD0_rz_simple, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_rz_simple" + sfx,
            .hists = {grouped(hists.real_pt2_MD1_rz_simple, grp, c), grouped(hists.fake_pt2_MD1_rz_simple, grp, c)},
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
            .hists = {grouped(hists.real_unused_pt2_deltaPT, grp, c), grouped(hists.fake_unused_pt2_deltaPT, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPHI" + sfx,
            .hists = {grouped(hists.real_unused_pt2_deltaPHI, grp, c), grouped(hists.fake_unused_pt2_deltaPHI, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Absolute #eta" + ttl,
            .xAxis = "#eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaETA" + sfx,
            .hists = {grouped(hists.real_unused_pt2_deltaETA, grp, c), grouped(hists.fake_unused_pt2_deltaETA, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Delta R" + ttl,
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaR" + sfx,
            .hists = {grouped(hists.real_unused_pt2_deltaR, grp, c), grouped(hists.fake_unused_pt2_deltaR, grp, c)},
            .legend = {"Real", "Fake"},
           .printYields = true
        });


        recipes.push_back({
            .title = "Unused LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdBeta" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTdBeta, grp, c), grouped(hists.fake_unused_pt2_LSTdBeta, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTKinZRes" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTKinZRes, grp, c), grouped(hists.fake_unused_pt2_LSTKinZRes, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTOrgZRes" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTOrgZRes, grp, c), grouped(hists.fake_unused_pt2_LSTOrgZRes, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST #Delta#phi" + ttl,
            .xAxis = "LST #Delta#phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdPhi" + sfx,
            .hists = {grouped(hists.real_unused_pt2_LSTdPhi, grp, c), grouped(hists.fake_unused_pt2_LSTdPhi, grp, c)},
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
            .hists = {grouped(hists.real_unused_pt2_MD0_dXY, grp, c), grouped(hists.fake_unused_pt2_MD0_dXY, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dz" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD0_dZ, grp, c), grouped(hists.fake_unused_pt2_MD0_dZ, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dxy" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD1_dXY, grp, c), grouped(hists.fake_unused_pt2_MD1_dXY, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dz" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD1_dZ, grp, c), grouped(hists.fake_unused_pt2_MD1_dZ, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_rz_simple" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD0_rz_simple, grp, c), grouped(hists.fake_unused_pt2_MD0_rz_simple, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_rz_simple" + sfx,
            .hists = {grouped(hists.real_unused_pt2_MD1_rz_simple, grp, c), grouped(hists.fake_unused_pt2_MD1_rz_simple, grp, c)},
            .legend = {"Real", "Fake"},
            .printYields = true
        });
     }
    }

    return recipes;
}
