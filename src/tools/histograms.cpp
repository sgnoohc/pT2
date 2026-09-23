#include "histograms.h"
#include <string>

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
}
