#include "extra_cuts.h"
#include <TMath.h>
#include <cmath>

namespace extra_cuts{
     
   /*  int getCategoryFromDetId(uint32_t detId) {
        int subdet = (detId >> 25) & 0x7;
        
        if (subdet == 4) { // ENDCAP
            int disk = (detId >> 18) & 0x7; 
            int ring = (detId >> 12) & 0xF;  // Extract the Ring number!
            bool isPS = (ring <= 10);        // Rings 1-10 are PS, 11-15 are 2S
            
            if (disk == 1) return isPS ? 7 : 8;
            if (disk == 2) return isPS ? 9 : 10;
            return isPS ? 11 : 12;           // Disk 3+ PS and 2S
        }
        if (subdet == 5) { // BARREL
            int layer = (detId >> 20) & 0x7; 
            int side  = (detId >> 18) & 0x3; 
            bool isTilted = (side == 1 || side == 2);
            
            if (layer == 1) return isTilted ? 1 : 0; // 0: L1 Flat, 1: L1 Tilted
            if (layer == 2) return isTilted ? 3 : 2; // 2: L2 Flat, 3: L2 Tilted
            if (layer == 3) return isTilted ? 5 : 4; // 4: L3 Flat, 5: L3 Tilted
            return 6; // Layers 4, 5, 6 (These are all 2S Flat)
        }
        return 6; // Fallback
    }*/
    int getCategoryFromDetId(uint32_t detId) {
        int subdet = (detId >> 25) & 0x7;

        if (subdet == 4) { // ENDCAP
            int disk = (detId >> 18) & 0x7;
            int ring = (detId >> 12) & 0xF;
            bool isPS = (ring <= 10);

            if (disk == 1) return isPS ? 7 : 8;
            if (disk == 2) return isPS ? 9 : 10;
            if (disk == 3) return isPS ? 11 : 12;  // STRICTLY Disk 3 ONLY

            return -1; // Disks 4 and 5 are completely rejected
        }
        if (subdet == 5) { // BARREL
            int layer = (detId >> 20) & 0x7;
            int side  = (detId >> 18) & 0x3;
            bool isTilted = (side == 1 || side == 2);

            if (layer == 1) return isTilted ? 1 : 0;
            if (layer == 2) return isTilted ? 3 : 2;
            if (layer == 3) return isTilted ? 5 : 4;

            return -1; // Layers 4, 5, and 6 are completely rejected
        }
        return -1; // Fallback for safety
    }

    int getConnectionIndex(int c0, int c1) {
        if (c0 == 0 && c1 == 2) return 0;   // Bar_L1_F -> Bar_L2_F
        if (c0 == 0 && c1 == 3) return 1;   // Bar_L1_F -> Bar_L2_T
        if (c0 == 1 && c1 == 2) return 2;   // Bar_L1_T -> Bar_L2_F
        if (c0 == 1 && c1 == 3) return 3;   // Bar_L1_T -> Bar_L2_T
        if (c0 == 1 && c1 == 7) return 4;   // Bar_L1_T -> Enc_D1_PS

        if (c0 == 2 && c1 == 4) return 5;   // Bar_L2_F -> Bar_L3_F
        if (c0 == 2 && c1 == 5) return 6;   // Bar_L2_F -> Bar_L3_T
        if (c0 == 3 && c1 == 4) return 7;   // Bar_L2_T -> Bar_L3_F
        if (c0 == 3 && c1 == 5) return 8;   // Bar_L2_T -> Bar_L3_T
        if (c0 == 3 && c1 == 7) return 9;   // Bar_L2_T -> Enc_D1_PS

        if (c0 == 7 && c1 == 9) return 10;  // Enc_D1_PS -> Enc_D2_PS
        if (c0 == 7 && c1 == 10) return 11; // Enc_D1_PS -> Enc_D2_2S
        if (c0 == 9 && c1 == 11) return 12; // Enc_D2_PS -> Enc_D3p_PS
        return -1; // Invalid connection
    }
    
    // =========================================================================
    //                        ALPAKA MATH REPLICAS
    // =========================================================================

    // 1. Exact replica of cms::alpakatools::reducePhiRange
    double reducePhiRange(double x) {
        constexpr double o2pi = 1.0 / (2.0 * TMath::Pi());
        if (std::abs(x) <= TMath::Pi()) return x;
        double n = std::round(x * o2pi);
        return x - n * 2.0 * TMath::Pi();
    }

    // 2. Exact replica of cms::alpakatools::deltaPhi
    double deltaPhiLST(double x1, double y1, double x2, double y2) {
        return reducePhiRange(std::atan2(-y2, -x2) - std::atan2(-y1, -x1));
    }

    // =========================================================================
    //                        LST KINEMATIC VARIABLES
    // =========================================================================

    // 3. The actual dPhi calculation from runTripletDefaultAlgoPPBB
    double calculateLSTDPhi(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -900.0f; 

        // Get InLo (Inner hit of the Pixel Segment)
        double x_InLo = data.pls_hit0_x->at(pls_idx);
        double y_InLo = data.pls_hit0_y->at(pls_idx);

        // Fallback: If hit0 is missing, use hit1
        if (x_InLo < invalid_val) {
            x_InLo = data.pls_hit1_x->at(pls_idx);
            y_InLo = data.pls_hit1_y->at(pls_idx);
        }
        // If still invalid, return dummy value
        if (x_InLo < invalid_val) return -999.0; 

        // Get OutLo (Inner MD of the Outer Tracker Segment)
        int md0_idx = data.ls_mdIdx0->at(ls_idx);
        if (md0_idx < 0) return -999.0;
        
        double x_OutLo = data.md_anchor_x->at(md0_idx);
        double y_OutLo = data.md_anchor_y->at(md0_idx);

        // Replicate the exact LST algorithm logic:
        double midPointX = 0.5 * (x_InLo + x_OutLo);
        double midPointY = 0.5 * (y_InLo + y_OutLo);

        double diffX = x_OutLo - x_InLo;
        double diffY = y_OutLo - y_InLo;

        return deltaPhiLST(midPointX, midPointY, diffX, diffY);
    }
    std::vector<double> calculateLSTdBeta(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -900.0f;

        // If the hits are invalid, return dummy values
        std::vector<double> invalid_result = {-999.0, -999.0, -999.0};

        const double k2Rinv1GeVf = 0.00570361;
        const double kSinAlphaMax = 0.95;

        double x_InUp = data.pls_hit3_x->at(pls_idx);
        double y_InUp = data.pls_hit3_y->at(pls_idx);

        if (x_InUp < invalid_val) {
            x_InUp = data.pls_hit2_x->at(pls_idx);
            y_InUp = data.pls_hit2_y->at(pls_idx);
        }
        if (x_InUp < invalid_val) return invalid_result;

        double pt = data.pls_pt->at(pls_idx);
        double phi = data.pls_phi->at(pls_idx);
        double px = pt * std::cos(phi);
        double py = pt * std::sin(phi);

        int md0_idx = data.ls_mdIdx0->at(ls_idx);
        int md1_idx = data.ls_mdIdx1->at(ls_idx);
        if (md0_idx < 0 || md1_idx < 0) return invalid_result;

        double x_OutLo = data.md_anchor_x->at(md0_idx);
        double y_OutLo = data.md_anchor_y->at(md0_idx);
        double x_OutUp = data.md_anchor_x->at(md1_idx);
        double y_OutUp = data.md_anchor_y->at(md1_idx);

        double tl_axis_x = x_OutUp - x_InUp;
        double tl_axis_y = y_OutUp - y_InUp;

        double betaIn = -deltaPhiLST(px, py, tl_axis_x, tl_axis_y);

        double alpha_OutUp = deltaPhiLST(x_OutUp, y_OutUp, x_OutUp - x_OutLo, y_OutUp - y_OutLo);
        double betaOut = -alpha_OutUp + deltaPhiLST(x_OutUp, y_OutUp, tl_axis_x, tl_axis_y);

        // LST Curved Path Correction (for lIn = 0)
        double sdOut_dr = std::sqrt(std::pow(x_OutUp - x_OutLo, 2) + std::pow(y_OutUp - y_OutLo, 2));
        double abs_pt = std::abs(pt);

        double correction_arg = std::min((sdOut_dr * k2Rinv1GeVf) / abs_pt, kSinAlphaMax);
        double correction = std::asin(correction_arg);
        betaOut += std::copysign(correction, betaOut);

        double dBeta = betaIn - betaOut;

        // Return them as a simple array/vector!
        // Index 0: betaIn
        // Index 1: betaOut
        // Index 2: dBeta
        return {betaIn, betaOut, dBeta};
    }
    // =========================================================================
    //                  Z-AXIS WINDOWS (ORIGIN & KINEMATIC)
    // =========================================================================

    // 1. KINEMATIC Z-RESIDUAL (zLoPointed / zHiPointed)
    double calculateLSTKinematicZResidual(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -900.0f;
        const double k2Rinv1GeVf = 0.00570361;

        double x_InUp = data.pls_hit3_x->at(pls_idx);
        double y_InUp = data.pls_hit3_y->at(pls_idx);
        double z_InUp = data.pls_hit3_z->at(pls_idx);

        if (x_InUp < invalid_val) {
            x_InUp = data.pls_hit2_x->at(pls_idx);
            y_InUp = data.pls_hit2_y->at(pls_idx);
            z_InUp = data.pls_hit2_z->at(pls_idx);
        }
        if (x_InUp < invalid_val) return -999.0;

        double rt_InUp = std::sqrt(x_InUp*x_InUp + y_InUp*y_InUp);

        int md0_idx = data.ls_mdIdx0->at(ls_idx);
        if (md0_idx < 0) return -999.0;

        double x_OutLo = data.md_anchor_x->at(md0_idx);
        double y_OutLo = data.md_anchor_y->at(md0_idx);
        double z_OutLo = data.md_anchor_z->at(md0_idx);
        double rt_OutLo = std::sqrt(x_OutLo*x_OutLo + y_OutLo*y_OutLo);

        double ptIn = data.pls_pt->at(pls_idx);
        double eta = data.pls_eta->at(pls_idx);

        // pz / ptIn = sinh(eta)
        double dzDrIn = std::sinh(eta);

        double drt_OutLo_InUp = rt_OutLo - rt_InUp;
        double curve_corr = (drt_OutLo_InUp * drt_OutLo_InUp * 4.0 * k2Rinv1GeVf * k2Rinv1GeVf) / (ptIn * ptIn * 24.0);

        double dzMean = dzDrIn * drt_OutLo_InUp * (1.0 + curve_corr);
        double z_OutLo_predicted = z_InUp + dzMean;

        return (z_OutLo - z_OutLo_predicted);
    }

    // 2. ORIGIN GEOMETRIC Z-RESIDUAL (zLo / zHi)
    double calculateLSTOriginZResidual(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -900.0f;

        double x_InUp = data.pls_hit3_x->at(pls_idx);
        double y_InUp = data.pls_hit3_y->at(pls_idx);
        double z_InUp = data.pls_hit3_z->at(pls_idx);

        if (x_InUp < invalid_val) {
            x_InUp = data.pls_hit2_x->at(pls_idx);
            y_InUp = data.pls_hit2_y->at(pls_idx);
            z_InUp = data.pls_hit2_z->at(pls_idx);
        }
        if (x_InUp < invalid_val) return -999.0;

        double rt_InUp = std::sqrt(x_InUp*x_InUp + y_InUp*y_InUp);

        int md0_idx = data.ls_mdIdx0->at(ls_idx);
        if (md0_idx < 0) return -999.0;

        double x_OutLo = data.md_anchor_x->at(md0_idx);
        double y_OutLo = data.md_anchor_y->at(md0_idx);
        double z_OutLo = data.md_anchor_z->at(md0_idx);
        double rt_OutLo = std::sqrt(x_OutLo*x_OutLo + y_OutLo*y_OutLo);

        // Simple straight-line projection from Origin (0,0,0)
        double z_OutLo_predicted = (z_InUp / rt_InUp) * rt_OutLo;

        return (z_OutLo - z_OutLo_predicted);
    }

} 
