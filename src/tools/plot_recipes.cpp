#include "plot_recipes.h"
#include "histograms.h"
#include <vector>
#include <string>

std::vector<PlotRecipe> getPt2Recipes(const HistogramManager& hists) {
    std::vector<PlotRecipe> recipes;

    for (int i = 0; i < 13; ++i) {
        for (int c = 0; c < 2; ++c) {
        std::string sfx = "_" + hists.catNames[i] + "_" + hists.chargeNames[c];
        std::string ttl = " (" + hists.catTitles[i] +  ", " + hists.chargeTitles[c] + ")";

        // =====================================================================
        // ALL pT2s - KINEMATICS & LST VARIABLES
        // =====================================================================

        recipes.push_back({
            .title = "pT2 Delta p_{T}" + ttl,
            .xAxis = "#Delta p_{T} [GeV]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPT" + sfx,
            .hists = {hists.real_pt2_deltaPT[i][c], hists.fake_pt2_deltaPT[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaPHI" + sfx,
            .hists = {hists.real_pt2_deltaPHI[i][c], hists.fake_pt2_deltaPHI[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta #eta" + ttl,
            .xAxis = "#eta",
            .yAxis = "Entries",
            .filename = "pt2_all_ETA" + sfx,
            .hists = {hists.real_pt2_deltaETA[i][c], hists.fake_pt2_deltaETA[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // --- pLS Absolute Eta ---
        recipes.push_back({
            .title = "pLS Absolute #eta" + ttl,
            .xAxis = "pLS #eta",
            .yAxis = "Entries",
            .filename = "pt2_all_pls_ETA" + sfx,
            .hists = {hists.real_pt2_pls_ETA[i][c], hists.fake_pt2_pls_ETA[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pLS Absolute #eta" + ttl,
            .xAxis = "pLS #eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_pls_ETA" + sfx,
            .hists = {hists.real_unused_pt2_pls_ETA[i][c], hists.fake_unused_pt2_pls_ETA[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        // --- LS Absolute Eta ---
        recipes.push_back({
            .title = "LS Absolute #eta" + ttl,
            .xAxis = "LS #eta",
            .yAxis = "Entries",
            .filename = "pt2_all_ls_ETA" + sfx,
            .hists = {hists.real_pt2_ls_ETA[i][c], hists.fake_pt2_ls_ETA[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LS Absolute #eta" + ttl,
            .xAxis = "LS #eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_ls_ETA" + sfx,
            .hists = {hists.real_unused_pt2_ls_ETA[i][c], hists.fake_unused_pt2_ls_ETA[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "pT2 Delta R" + ttl,
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_all_deltaR" + sfx,
            .hists = {hists.real_pt2_deltaR[i][c], hists.fake_pt2_deltaR[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdBeta" + sfx,
            .hists = {hists.real_pt2_LSTdBeta[i][c], hists.fake_pt2_LSTdBeta[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTKinZRes" + sfx,
            .hists = {hists.real_pt2_LSTKinZRes[i][c], hists.fake_pt2_LSTKinZRes[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTOrgZRes" + sfx,
            .hists = {hists.real_pt2_LSTOrgZRes[i][c], hists.fake_pt2_LSTOrgZRes[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "LST #Delta#phi" + ttl,
            .xAxis = "LST #Delta#phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_all_LSTdPhi" + sfx,
            .hists = {hists.real_pt2_LSTdPhi[i][c], hists.fake_pt2_LSTdPhi[i][c]},
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
            .hists = {hists.real_pt2_MD0_dXY[i][c], hists.fake_pt2_MD0_dXY[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_dz" + sfx,
            .hists = {hists.real_pt2_MD0_dZ[i][c], hists.fake_pt2_MD0_dZ[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dxy" + sfx,
            .hists = {hists.real_pt2_MD1_dXY[i][c], hists.fake_pt2_MD1_dXY[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_dz" + sfx,
            .hists = {hists.real_pt2_MD1_dZ[i][c], hists.fake_pt2_MD1_dZ[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md0_rz_simple" + sfx,
            .hists = {hists.real_pt2_MD0_rz_simple[i][c], hists.fake_pt2_MD0_rz_simple[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_all_md1_rz_simple" + sfx,
            .hists = {hists.real_pt2_MD1_rz_simple[i][c], hists.fake_pt2_MD1_rz_simple[i][c]},
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
            .hists = {hists.real_unused_pt2_deltaPT[i][c], hists.fake_unused_pt2_deltaPT[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Delta #phi" + ttl,
            .xAxis = "#Delta #phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaPHI" + sfx,
            .hists = {hists.real_unused_pt2_deltaPHI[i][c], hists.fake_unused_pt2_deltaPHI[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Absolute #eta" + ttl,
            .xAxis = "#eta",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaETA" + sfx,
            .hists = {hists.real_unused_pt2_deltaETA[i][c], hists.fake_unused_pt2_deltaETA[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused pT2 Delta R" + ttl,
            .xAxis = "#Delta R",
            .yAxis = "Entries",
            .filename = "pt2_unused_deltaR" + sfx,
            .hists = {hists.real_unused_pt2_deltaR[i][c], hists.fake_unused_pt2_deltaR[i][c]},
            .legend = {"Real", "Fake"},
           .printYields = true
        });


        recipes.push_back({
            .title = "Unused LST Delta Beta" + ttl,
            .xAxis = "LST #Delta#beta [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdBeta" + sfx,
            .hists = {hists.real_unused_pt2_LSTdBeta[i][c], hists.fake_unused_pt2_LSTdBeta[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Kinematic Z-Residual" + ttl,
            .xAxis = "Actual Z - Kinematic Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTKinZRes" + sfx,
            .hists = {hists.real_unused_pt2_LSTKinZRes[i][c], hists.fake_unused_pt2_LSTKinZRes[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST Geometric Z-Residual" + ttl,
            .xAxis = "Actual Z - Origin Predicted Z [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTOrgZRes" + sfx,
            .hists = {hists.real_unused_pt2_LSTOrgZRes[i][c], hists.fake_unused_pt2_LSTOrgZRes[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused LST #Delta#phi" + ttl,
            .xAxis = "LST #Delta#phi [rad]",
            .yAxis = "Entries",
            .filename = "pt2_unused_LSTdPhi" + sfx,
            .hists = {hists.real_unused_pt2_LSTdPhi[i][c], hists.fake_unused_pt2_LSTdPhi[i][c]},
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
            .hists = {hists.real_unused_pt2_MD0_dXY[i][c], hists.fake_unused_pt2_MD0_dXY[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_dz" + sfx,
            .hists = {hists.real_unused_pt2_MD0_dZ[i][c], hists.fake_unused_pt2_MD0_dZ[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Transverse Distance (dXY)" + ttl,
            .xAxis = "dXY [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dxy" + sfx,
            .hists = {hists.real_unused_pt2_MD1_dXY[i][c], hists.fake_unused_pt2_MD1_dXY[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 Longitudinal Distance (dZ)" + ttl,
            .xAxis = "dZ [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_dz" + sfx,
            .hists = {hists.real_unused_pt2_MD1_dZ[i][c], hists.fake_unused_pt2_MD1_dZ[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD0 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md0_rz_simple" + sfx,
            .hists = {hists.real_unused_pt2_MD0_rz_simple[i][c], hists.fake_unused_pt2_MD0_rz_simple[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });

        recipes.push_back({
            .title = "Unused MD1 R-Z Simple Residual" + ttl,
            .xAxis = "R_{act} - R_{pred} [cm]",
            .yAxis = "Entries",
            .filename = "pt2_unused_md1_rz_simple" + sfx,
            .hists = {hists.real_unused_pt2_MD1_rz_simple[i][c], hists.fake_unused_pt2_MD1_rz_simple[i][c]},
            .legend = {"Real", "Fake"},
            .printYields = true
        });
     }
    }

    return recipes;
}
