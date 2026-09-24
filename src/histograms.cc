#include "histograms.h"
#include <TFile.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>

void HistogramManager::init() {
    for (int i = 0; i < 13; ++i) {
        for (int c = 0; c < 2; ++c) {
        std::string sfx = "_" + catNames[i] + "_" + chargeNames[c];
        std::string ttl = " (" + catTitles[i] +  ", " + chargeTitles[c] + ")";

        // =====================================================================
        // KINEMATICS
        // =====================================================================

        real_pt2_deltaPT[i][c] = new TH1D(
            ("real_pt2_deltaPT" + sfx).c_str(),
            ("Real pT2 #Delta pT" + ttl + "; #Delta pT [GeV]; Entries").c_str(),
            180, -10.0, 10.0
        );

        fake_pt2_deltaPT[i][c] = new TH1D(
            ("fake_pt2_deltaPT" + sfx).c_str(),
            ("Fake pT2 #Delta pT" + ttl + "; #Delta pT [GeV]; Entries").c_str(),
            180, -10.0, 10.0
        );

        real_unused_pt2_deltaPT[i][c] = new TH1D(
            ("real_unused_pt2_deltaPT" + sfx).c_str(),
            ("Real Unused pT2 #Delta pT" + ttl + "; #Delta pT [GeV]; Entries").c_str(),
            180, -10.0, 10.0
        );

        fake_unused_pt2_deltaPT[i][c] = new TH1D(
            ("fake_unused_pt2_deltaPT" + sfx).c_str(),
            ("Fake Unused pT2 #Delta pT" + ttl + "; #Delta pT [GeV]; Entries").c_str(),
            180, -10.0, 10.0
        );

        real_pt2_deltaETA[i][c] = new TH1D(
            ("real_pt2_deltaETA" + sfx).c_str(),
            ("Real pT2 #Delta #eta" + ttl + "; #Delta #eta; Entries").c_str(),
            180, -1.0, 1.0
        );

        fake_pt2_deltaETA[i][c] = new TH1D(
            ("fake_pt2_deltaETA" + sfx).c_str(),
            ("Fake pT2 #Delta #eta" + ttl + "; #Delta #eta; Entries").c_str(),
            180, -1.0, 1.0
        );

        real_unused_pt2_deltaETA[i][c] = new TH1D(
            ("real_unused_pt2_deltaETA" + sfx).c_str(),
            ("Real Unused pT2 #Delta #eta" + ttl + "; #Delta #eta; Entries").c_str(),
            180, -1.0, 1.0
        );

        fake_unused_pt2_deltaETA[i][c] = new TH1D(
            ("fake_unused_pt2_deltaETA" + sfx).c_str(),
            ("Fake Unused pT2 #Delta #eta" + ttl + "; #Delta #eta; Entries").c_str(),
            180, -1.0, 1.0
        );

        real_pt2_deltaPHI[i][c] = new TH1D(
            ("real_pt2_deltaPHI" + sfx).c_str(),
            ("Real pT2 #Delta #phi" + ttl + "; #Delta #phi [rad]; Entries").c_str(),
            180, -1.0, 1.0
        );

        fake_pt2_deltaPHI[i][c] = new TH1D(
            ("fake_pt2_deltaPHI" + sfx).c_str(),
            ("Fake pT2 #Delta #phi" + ttl + "; #Delta #phi [rad]; Entries").c_str(),
            180, -1.0, 1.0
        );

        real_unused_pt2_deltaPHI[i][c] = new TH1D(
            ("real_unused_pt2_deltaPHI" + sfx).c_str(),
            ("Real Unused pT2 #Delta #phi" + ttl + "; #Delta #phi [rad]; Entries").c_str(),
            180, -1.0, 1.0
        );

        fake_unused_pt2_deltaPHI[i][c] = new TH1D(
            ("fake_unused_pt2_deltaPHI" + sfx).c_str(),
            ("Fake Unused pT2 #Delta #phi" + ttl + "; #Delta #phi [rad]; Entries").c_str(),
            180, -1.0, 1.0
        );

        // pLS Absolute Eta
        real_pt2_pls_ETA[i][c] = new TH1D(
            ("real_pt2_pls_ETA"+sfx).c_str(), 
            ("Real pLS #eta"+ttl+";pLS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );
        
        fake_pt2_pls_ETA[i][c] = new TH1D(
            ("fake_pt2_pls_ETA"+sfx).c_str(), 
            ("Fake pLS #eta"+ttl+";pLS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );
        
        real_unused_pt2_pls_ETA[i][c] = new TH1D(
            ("real_unused_pt2_pls_ETA"+sfx).c_str(), 
            ("Real Unused pLS #eta"+ttl+";pLS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );
        
        fake_unused_pt2_pls_ETA[i][c] = new TH1D(
            ("fake_unused_pt2_pls_ETA"+sfx).c_str(), 
            ("Fake Unused pLS #eta"+ttl+";pLS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );

        // LS Absolute Eta
        real_pt2_ls_ETA[i][c] = new TH1D(
            ("real_pt2_ls_ETA"+sfx).c_str(), 
            ("Real LS #eta"+ttl+";LS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );
        
        fake_pt2_ls_ETA[i][c] = new TH1D(
            ("fake_pt2_ls_ETA"+sfx).c_str(), 
            ("Fake LS #eta"+ttl+";LS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );
        
        real_unused_pt2_ls_ETA[i][c] = new TH1D(
            ("real_unused_pt2_ls_ETA"+sfx).c_str(), 
            ("Real Unused LS #eta"+ttl+";LS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );
        
        fake_unused_pt2_ls_ETA[i][c] = new TH1D(
            ("fake_unused_pt2_ls_ETA"+sfx).c_str(), 
            ("Fake Unused LS #eta"+ttl+";LS #eta;Entries").c_str(), 
            180, -4.0, 4.0
        );
        //-------------------------------------------------
        real_pt2_deltaR[i][c] = new TH1D(
            ("real_pt2_deltaR" + sfx).c_str(),
            ("Real pT2 #Delta R" + ttl + "; #Delta R; Entries").c_str(),
            180, 0, 1.0
        );

        fake_pt2_deltaR[i][c] = new TH1D(
            ("fake_pt2_deltaR" + sfx).c_str(),
            ("Fake pT2 #Delta R" + ttl + "; #Delta R; Entries").c_str(),
            180, 0, 1.0
        );

        real_unused_pt2_deltaR[i][c] = new TH1D(
            ("real_unused_pt2_deltaR" + sfx).c_str(),
            ("Real Unused pT2 #Delta R" + ttl + "; #Delta R; Entries").c_str(),
            180, 0, 1.0
        );

        fake_unused_pt2_deltaR[i][c] = new TH1D(
            ("fake_unused_pt2_deltaR" + sfx).c_str(),
            ("Fake Unused pT2 #Delta R" + ttl + "; #Delta R; Entries").c_str(),
            180, 0, 1.0
        );

        real_pt2_deltaAngle[i][c] = new TH1D(
            ("real_pt2_deltaAngle" + sfx).c_str(),
            ("Real pT2 Directional #Delta#alpha" + ttl + "; #Delta#alpha [rad]; Entries").c_str(),
            180, 0, 0.5
        );

        fake_pt2_deltaAngle[i][c] = new TH1D(
            ("fake_pt2_deltaAngle" + sfx).c_str(),
            ("Fake pT2 Directional #Delta#alpha" + ttl + "; #Delta#alpha [rad]; Entries").c_str(),
            180, 0, 0.5
        );

        real_unused_pt2_deltaAngle[i][c] = new TH1D(
            ("real_unused_pt2_deltaAngle" + sfx).c_str(),
            ("Real Unused pT2 #Delta#alpha" + ttl + "; #Delta#alpha [rad]; Entries").c_str(),
            180, 0, 0.5
        );

        fake_unused_pt2_deltaAngle[i][c] = new TH1D(
            ("fake_unused_pt2_deltaAngle" + sfx).c_str(),
            ("Fake Unused pT2 #Delta#alpha" + ttl + "; #Delta#alpha [rad]; Entries").c_str(),
            180, 0, 0.5
        );

        // =====================================================================
        // LST VARIABLES
        // =====================================================================

        real_pt2_LSTdPhi[i][c] = new TH1D(
            ("real_pt2_LSTdPhi" + sfx).c_str(),
            ("Real pT2 LST #Delta#phi" + ttl + "; LST #Delta#phi [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        fake_pt2_LSTdPhi[i][c] = new TH1D(
            ("fake_pt2_LSTdPhi" + sfx).c_str(),
            ("Fake pT2 LST #Delta#phi" + ttl + "; LST #Delta#phi [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        real_unused_pt2_LSTdPhi[i][c] = new TH1D(
            ("real_unused_pt2_LSTdPhi" + sfx).c_str(),
            ("Real Unused LST #Delta#phi" + ttl + "; LST #Delta#phi [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        fake_unused_pt2_LSTdPhi[i][c] = new TH1D(
            ("fake_unused_pt2_LSTdPhi" + sfx).c_str(),
            ("Fake Unused LST #Delta#phi" + ttl + "; LST #Delta#phi [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        real_pt2_LSTdBeta[i][c] = new TH1D(
            ("real_pt2_LSTdBeta" + sfx).c_str(),
            ("Real pT2 LST #Delta#beta" + ttl + "; LST #Delta#beta [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        fake_pt2_LSTdBeta[i][c] = new TH1D(
            ("fake_pt2_LSTdBeta" + sfx).c_str(),
            ("Fake pT2 LST #Delta#beta" + ttl + "; LST #Delta#beta [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        real_unused_pt2_LSTdBeta[i][c] = new TH1D(
            ("real_unused_pt2_LSTdBeta" + sfx).c_str(),
            ("Real Unused LST #Delta#beta" + ttl + "; LST #Delta#beta [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        fake_unused_pt2_LSTdBeta[i][c] = new TH1D(
            ("fake_unused_pt2_LSTdBeta" + sfx).c_str(),
            ("Fake Unused LST #Delta#beta" + ttl + "; LST #Delta#beta [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        real_pt2_LSTbetaOut[i][c] = new TH1D(
            ("real_pt2_LSTbetaOut" + sfx).c_str(),
            ("Real pT2 LST #beta_{Out}" + ttl + "; LST #beta_{Out} [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        fake_pt2_LSTbetaOut[i][c] = new TH1D(
            ("fake_pt2_LSTbetaOut" + sfx).c_str(),
            ("Fake pT2 LST #beta_{Out}" + ttl + "; LST #beta_{Out} [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        real_unused_pt2_LSTbetaOut[i][c] = new TH1D(
            ("real_unused_pt2_LSTbetaOut" + sfx).c_str(),
            ("Real Unused LST #beta_{Out}" + ttl + "; LST #beta_{Out} [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        fake_unused_pt2_LSTbetaOut[i][c] = new TH1D(
            ("fake_unused_pt2_LSTbetaOut" + sfx).c_str(),
            ("Fake Unused LST #beta_{Out}" + ttl + "; LST #beta_{Out} [rad]; Entries").c_str(),
            180, -0.5, 0.5
        );

        real_pt2_LSTKinZRes[i][c] = new TH1D(
            ("real_pt2_LSTKinZRes" + sfx).c_str(),
            ("Real pT2 Kinematic Z-Res" + ttl + "; Z_{act} - Z_{pred} [cm]; Entries").c_str(),
            180, -5.0, 5.0
        );

        fake_pt2_LSTKinZRes[i][c] = new TH1D(
            ("fake_pt2_LSTKinZRes" + sfx).c_str(),
            ("Fake pT2 Kinematic Z-Res" + ttl + "; Z_{act} - Z_{pred} [cm]; Entries").c_str(),
            180, -5.0, 5.0
        );

        real_unused_pt2_LSTKinZRes[i][c] = new TH1D(
            ("real_unused_pt2_LSTKinZRes" + sfx).c_str(),
            ("Real Unused Kinematic Z-Res" + ttl + "; Z_{act} - Z_{pred} [cm]; Entries").c_str(),
            180, -5.0, 5.0
        );

        fake_unused_pt2_LSTKinZRes[i][c] = new TH1D(
            ("fake_unused_pt2_LSTKinZRes" + sfx).c_str(),
            ("Fake Unused Kinematic Z-Res" + ttl + "; Z_{act} - Z_{pred} [cm]; Entries").c_str(),
            180, -5.0, 5.0
        );

        real_pt2_LSTOrgZRes[i][c] = new TH1D(
            ("real_pt2_LSTOrgZRes" + sfx).c_str(),
            ("Real pT2 Geometric Z-Res" + ttl + "; Z_{act} - Z_{org} [cm]; Entries").c_str(),
            180, -25.0, 25.0
        );

        fake_pt2_LSTOrgZRes[i][c] = new TH1D(
            ("fake_pt2_LSTOrgZRes" + sfx).c_str(),
            ("Fake pT2 Geometric Z-Res" + ttl + "; Z_{act} - Z_{org} [cm]; Entries").c_str(),
            180, -25.0, 25.0
        );

        real_unused_pt2_LSTOrgZRes[i][c] = new TH1D(
            ("real_unused_pt2_LSTOrgZRes" + sfx).c_str(),
            ("Real Unused Geometric Z-Res" + ttl + "; Z_{act} - Z_{org} [cm]; Entries").c_str(),
            180, -25.0, 25.0
        );

        fake_unused_pt2_LSTOrgZRes[i][c] = new TH1D(
            ("fake_unused_pt2_LSTOrgZRes" + sfx).c_str(),
            ("Fake Unused Geometric Z-Res" + ttl + "; Z_{act} - Z_{org} [cm]; Entries").c_str(),
            180, -25.0, 25.0
        );

        // =====================================================================
        // MD COMPONENTS (dXY, dZ, RZ)
        // =====================================================================

        real_pt2_MD0_dXY[i][c] = new TH1D(
            ("real_pt2_MD0_dXY" + sfx).c_str(),
            ("Real MD0 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        fake_pt2_MD0_dXY[i][c] = new TH1D(
            ("fake_pt2_MD0_dXY" + sfx).c_str(),
            ("Fake MD0 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        real_unused_pt2_MD0_dXY[i][c] = new TH1D(
            ("real_unused_pt2_MD0_dXY" + sfx).c_str(),
            ("Real Unused MD0 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        fake_unused_pt2_MD0_dXY[i][c] = new TH1D(
            ("fake_unused_pt2_MD0_dXY" + sfx).c_str(),
            ("Fake Unused MD0 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        real_pt2_MD0_dZ[i][c] = new TH1D(
            ("real_pt2_MD0_dZ" + sfx).c_str(),
            ("Real MD0 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        fake_pt2_MD0_dZ[i][c] = new TH1D(
            ("fake_pt2_MD0_dZ" + sfx).c_str(),
            ("Fake MD0 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        real_unused_pt2_MD0_dZ[i][c] = new TH1D(
            ("real_unused_pt2_MD0_dZ" + sfx).c_str(),
            ("Real Unused MD0 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        fake_unused_pt2_MD0_dZ[i][c] = new TH1D(
            ("fake_unused_pt2_MD0_dZ" + sfx).c_str(),
            ("Fake Unused MD0 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        real_pt2_MD1_dXY[i][c] = new TH1D(
            ("real_pt2_MD1_dXY" + sfx).c_str(),
            ("Real MD1 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        fake_pt2_MD1_dXY[i][c] = new TH1D(
            ("fake_pt2_MD1_dXY" + sfx).c_str(),
            ("Fake MD1 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        real_unused_pt2_MD1_dXY[i][c] = new TH1D(
            ("real_unused_pt2_MD1_dXY" + sfx).c_str(),
            ("Real Unused MD1 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        fake_unused_pt2_MD1_dXY[i][c] = new TH1D(
            ("fake_unused_pt2_MD1_dXY" + sfx).c_str(),
            ("Fake Unused MD1 #DeltaXY" + ttl + "; #DeltaXY [cm]; Entries").c_str(),
            180, 0, 5.0
        );

        real_pt2_MD1_dZ[i][c] = new TH1D(
            ("real_pt2_MD1_dZ" + sfx).c_str(),
            ("Real MD1 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        fake_pt2_MD1_dZ[i][c] = new TH1D(
            ("fake_pt2_MD1_dZ" + sfx).c_str(),
            ("Fake MD1 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        real_unused_pt2_MD1_dZ[i][c] = new TH1D(
            ("real_unused_pt2_MD1_dZ" + sfx).c_str(),
            ("Real Unused MD1 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        fake_unused_pt2_MD1_dZ[i][c] = new TH1D(
            ("fake_unused_pt2_MD1_dZ" + sfx).c_str(),
            ("Fake Unused MD1 #DeltaZ" + ttl + "; #DeltaZ [cm]; Entries").c_str(),
            180, 0, 10.0
        );

        real_pt2_MD0_rz_simple[i][c] = new TH1D(
            ("real_pt2_MD0_rz_simple" + sfx).c_str(),
            ("Real MD0 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );

        fake_pt2_MD0_rz_simple[i][c] = new TH1D(
            ("fake_pt2_MD0_rz_simple" + sfx).c_str(),
            ("Fake MD0 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );

        real_unused_pt2_MD0_rz_simple[i][c] = new TH1D(
            ("real_unused_pt2_MD0_rz_simple" + sfx).c_str(),
            ("Real Unused MD0 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );

        fake_unused_pt2_MD0_rz_simple[i][c] = new TH1D(
            ("fake_unused_pt2_MD0_rz_simple" + sfx).c_str(),
            ("Fake Unused MD0 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );

        real_pt2_MD1_rz_simple[i][c] = new TH1D(
            ("real_pt2_MD1_rz_simple" + sfx).c_str(),
            ("Real MD1 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );

        fake_pt2_MD1_rz_simple[i][c] = new TH1D(
            ("fake_pt2_MD1_rz_simple" + sfx).c_str(),
            ("Fake MD1 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );

        real_unused_pt2_MD1_rz_simple[i][c] = new TH1D(
            ("real_unused_pt2_MD1_rz_simple" + sfx).c_str(),
            ("Real Unused MD1 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );

        fake_unused_pt2_MD1_rz_simple[i][c] = new TH1D(
            ("fake_unused_pt2_MD1_rz_simple" + sfx).c_str(),
            ("Fake Unused MD1 R-Z Residual" + ttl + "; [cm]; Entries").c_str(),
            180, -10, 10
        );
    }
    }

    for (auto& [name, grid] : grids())
        for (int i = 0; i < kNCat; ++i)
            for (int c = 0; c < kNCharge; ++c)
                (*grid)[i][c]->SetDirectory(nullptr);

    buildSets();
}

std::vector<std::pair<std::string, HistGrid*>> HistogramManager::grids() {
#define REG(x) {#x, &x}
    return {
        REG(real_pt2_deltaPT),
        REG(fake_pt2_deltaPT),
        REG(real_unused_pt2_deltaPT),
        REG(fake_unused_pt2_deltaPT),
        REG(real_pt2_deltaETA),
        REG(fake_pt2_deltaETA),
        REG(real_unused_pt2_deltaETA),
        REG(fake_unused_pt2_deltaETA),
        REG(real_pt2_deltaPHI),
        REG(fake_pt2_deltaPHI),
        REG(real_unused_pt2_deltaPHI),
        REG(fake_unused_pt2_deltaPHI),
        REG(real_pt2_deltaR),
        REG(fake_pt2_deltaR),
        REG(real_unused_pt2_deltaR),
        REG(fake_unused_pt2_deltaR),
        REG(real_pt2_deltaAngle),
        REG(fake_pt2_deltaAngle),
        REG(real_unused_pt2_deltaAngle),
        REG(fake_unused_pt2_deltaAngle),
        REG(real_pt2_pls_ETA),
        REG(fake_pt2_pls_ETA),
        REG(real_unused_pt2_pls_ETA),
        REG(fake_unused_pt2_pls_ETA),
        REG(real_pt2_ls_ETA),
        REG(fake_pt2_ls_ETA),
        REG(real_unused_pt2_ls_ETA),
        REG(fake_unused_pt2_ls_ETA),
        REG(real_pt2_LSTdPhi),
        REG(fake_pt2_LSTdPhi),
        REG(real_unused_pt2_LSTdPhi),
        REG(fake_unused_pt2_LSTdPhi),
        REG(real_pt2_LSTdBeta),
        REG(fake_pt2_LSTdBeta),
        REG(real_unused_pt2_LSTdBeta),
        REG(fake_unused_pt2_LSTdBeta),
        REG(real_pt2_LSTbetaOut),
        REG(fake_pt2_LSTbetaOut),
        REG(real_unused_pt2_LSTbetaOut),
        REG(fake_unused_pt2_LSTbetaOut),
        REG(real_pt2_LSTOrgZRes),
        REG(fake_pt2_LSTOrgZRes),
        REG(real_unused_pt2_LSTOrgZRes),
        REG(fake_unused_pt2_LSTOrgZRes),
        REG(real_pt2_LSTKinZRes),
        REG(fake_pt2_LSTKinZRes),
        REG(real_unused_pt2_LSTKinZRes),
        REG(fake_unused_pt2_LSTKinZRes),
        REG(real_pt2_MD0_dXY),
        REG(fake_pt2_MD0_dXY),
        REG(real_unused_pt2_MD0_dXY),
        REG(fake_unused_pt2_MD0_dXY),
        REG(real_pt2_MD0_dZ),
        REG(fake_pt2_MD0_dZ),
        REG(real_unused_pt2_MD0_dZ),
        REG(fake_unused_pt2_MD0_dZ),
        REG(real_pt2_MD1_dXY),
        REG(fake_pt2_MD1_dXY),
        REG(real_unused_pt2_MD1_dXY),
        REG(fake_unused_pt2_MD1_dXY),
        REG(real_pt2_MD1_dZ),
        REG(fake_pt2_MD1_dZ),
        REG(real_unused_pt2_MD1_dZ),
        REG(fake_unused_pt2_MD1_dZ),
        REG(real_pt2_MD0_rz_simple),
        REG(fake_pt2_MD0_rz_simple),
        REG(real_unused_pt2_MD0_rz_simple),
        REG(fake_unused_pt2_MD0_rz_simple),
        REG(real_pt2_MD1_rz_simple),
        REG(fake_pt2_MD1_rz_simple),
        REG(real_unused_pt2_MD1_rz_simple),
        REG(fake_unused_pt2_MD1_rz_simple),
    };
#undef REG
}

void HistogramManager::buildSets() {
    std::unordered_map<std::string, HistGrid*> g;
    for (auto& [name, grid] : grids()) g[name] = grid;

    for (int real = 0; real < 2; ++real) {
        for (int unused = 0; unused < 2; ++unused) {
            std::string prefix = std::string(real ? "real" : "fake") + (unused ? "_unused" : "") + "_pt2_";
            for (int i = 0; i < kNCat; ++i) {
                for (int c = 0; c < kNCharge; ++c) {
                    Pt2HistSet& s = sets_[real ? 0 : 1][unused][i][c];
                    s.deltaPT = (*g[prefix + "deltaPT"])[i][c];
                    s.deltaETA = (*g[prefix + "deltaETA"])[i][c];
                    s.deltaPHI = (*g[prefix + "deltaPHI"])[i][c];
                    s.deltaR = (*g[prefix + "deltaR"])[i][c];
                    s.deltaAngle = (*g[prefix + "deltaAngle"])[i][c];
                    s.pls_ETA = (*g[prefix + "pls_ETA"])[i][c];
                    s.ls_ETA = (*g[prefix + "ls_ETA"])[i][c];
                    s.LSTdPhi = (*g[prefix + "LSTdPhi"])[i][c];
                    s.LSTdBeta = (*g[prefix + "LSTdBeta"])[i][c];
                    s.LSTbetaOut = (*g[prefix + "LSTbetaOut"])[i][c];
                    s.LSTOrgZRes = (*g[prefix + "LSTOrgZRes"])[i][c];
                    s.LSTKinZRes = (*g[prefix + "LSTKinZRes"])[i][c];
                    s.MD0_dXY = (*g[prefix + "MD0_dXY"])[i][c];
                    s.MD0_dZ = (*g[prefix + "MD0_dZ"])[i][c];
                    s.MD1_dXY = (*g[prefix + "MD1_dXY"])[i][c];
                    s.MD1_dZ = (*g[prefix + "MD1_dZ"])[i][c];
                    s.MD0_rz_simple = (*g[prefix + "MD0_rz_simple"])[i][c];
                    s.MD1_rz_simple = (*g[prefix + "MD1_rz_simple"])[i][c];
                }
            }
        }
    }
}

void HistogramManager::write(const std::string& path) {
    TFile f(path.c_str(), "RECREATE");
    if (f.IsZombie()) throw std::runtime_error("cannot create " + path);
    for (auto& [name, grid] : grids())
        for (int i = 0; i < kNCat; ++i)
            for (int c = 0; c < kNCharge; ++c)
                f.WriteObject((*grid)[i][c], (*grid)[i][c]->GetName());
    f.Close();
}

void HistogramManager::load(const std::string& path) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) throw std::runtime_error("cannot open " + path);
    for (auto& [name, grid] : grids()) {
        for (int i = 0; i < kNCat; ++i) {
            for (int c = 0; c < kNCharge; ++c) {
                std::string key = name + "_" + catNames[i] + "_" + chargeNames[c];
                TH1D* h = f->Get<TH1D>(key.c_str());
                if (!h) throw std::runtime_error("histogram " + key + " missing in " + path);
                h->SetDirectory(nullptr);
                (*grid)[i][c] = h;
            }
        }
    }
    buildSets();
}
