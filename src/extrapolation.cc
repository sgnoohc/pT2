#include "extrapolation.h"
#include "rootReader.h"
// ROOT Includes
#include <TVector3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/Functor.h>
#include <cmath>
#include <iostream>

namespace extrapolation {

    // =========================================================================
    //                            HELPER FUNCTIONS
    // =========================================================================

    ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double R) {
        if (hits.size() < 2) return {0, 0};

        ROOT::Math::XYVector p1 = hits.front();
        ROOT::Math::XYVector p2 = hits.back();

        double dx = p2.X() - p1.X();
        double dy = p2.Y() - p1.Y();
        double d2 = dx*dx + dy*dy;
        double d = std::sqrt(d2);

        if (d > 2.0 * R || d == 0) return {0, 0};

        double midX = (p1.X() + p2.X()) / 2.0;
        double midY = (p1.Y() + p2.Y()) / 2.0;

        double h = std::sqrt(std::max(0.0, R*R - d2/4.0));

        double c1x = midX + h * (-dy / d);
        double c1y = midY + h * (dx / d);
        double c2x = midX - h * (-dy / d);
        double c2y = midY - h * (dx / d);

        if (hits.size() > 2) {
            ROOT::Math::XYVector pMid = hits[hits.size()/2];
            auto distSq = [](double x1, double y1, double x2, double y2) {
                return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
            };
            double diff1 = std::abs(std::sqrt(distSq(pMid.X(), pMid.Y(), c1x, c1y)) - R);
            double diff2 = std::abs(std::sqrt(distSq(pMid.X(), pMid.Y(), c2x, c2y)) - R);
            return (diff1 < diff2) ? ROOT::Math::XYVector(c1x, c1y) : ROOT::Math::XYVector(c2x, c2y);
        }

        return {c1x, c1y};
    }

std::vector<double> extrapolatePlsHelicallyAndGetDistance(int pls_idx, int ls_idx, const rootReader& data) {
    // Return format: {dXY_MD0, dZ_MD0, dXY_MD1, dZ_MD1}
    const std::vector<double> invalid_return = {-1.0, -999.0, -1.0, -999.0};
    const float invalid_hit_val = -900.0f;

    // 1. Collect valid pLS hits to find the helical pitch (Z vs Phi)
    std::vector<ROOT::Math::XYVector> hits_xy;
    std::vector<double> hits_z;
    auto add_hit = [&](float x, float y, float z) {
        if (x > invalid_hit_val) {
            hits_xy.emplace_back(x, y);
            hits_z.push_back(z);
        }
    };

    add_hit(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx));
    add_hit(data.pls_hit1_x->at(pls_idx), data.pls_hit1_y->at(pls_idx), data.pls_hit1_z->at(pls_idx));
    add_hit(data.pls_hit2_x->at(pls_idx), data.pls_hit2_y->at(pls_idx), data.pls_hit2_z->at(pls_idx));
    add_hit(data.pls_hit3_x->at(pls_idx), data.pls_hit3_y->at(pls_idx), data.pls_hit3_z->at(pls_idx));

    if (hits_xy.size() < 2) return invalid_return;

    // 2. Use official Circle Center and calculate Radius from pT
 //   double centerX = data.pLS_circleCenterX->at(pls_idx);
//    double centerY = data.pLS_circleCenterY->at(pls_idx);
//    ROOT::Math::XYVector circle_center(centerX, centerY);
   
    float reco_pt = data.pls_pt->at(pls_idx);
    double R_pred_cm = (reco_pt / (0.3 * 3.8)) * 100.0;
    
    ROOT::Math::XYVector circle_center = fitCircleWithFixedRadius(hits_xy, R_pred_cm);
    if (circle_center.X() == 0 && circle_center.Y() == 0) return invalid_return;

    // 3. Simple Linear Z-Pitch (z = a*phi + b)
    double first_angle = 0, first_z = 0;
    double last_angle = 0, last_z = 0;

    for (size_t i = 0; i < hits_xy.size(); ++i) {
        double angle = TMath::ATan2(hits_xy[i].Y() - circle_center.Y(), hits_xy[i].X() - circle_center.X());
        if (i == 0) {
            first_angle = angle;
            first_z = hits_z[i];
            last_angle = angle;
        } else {
            // Unwrap angle to ensure the helix doesn't "jump" 2pi
            while (angle - last_angle > TMath::Pi()) angle -= 2 * TMath::Pi();
            while (angle - last_angle < -TMath::Pi()) angle += 2 * TMath::Pi();
            last_angle = angle;
            last_z = hits_z[i];
        }
    }

    double delta_phi_pLS = last_angle - first_angle;
    if (std::abs(delta_phi_pLS) < 1e-9) return invalid_return;

    double a = (last_z - first_z) / delta_phi_pLS; // dz/dphi
    double b = first_z - a * first_angle;          // intercept

    // 4. Target Mini-Doublets (LS)
    int md0 = data.ls_mdIdx0->at(ls_idx);
    int md1 = data.ls_mdIdx1->at(ls_idx);
    if (md0 < 0 || md1 < 0) return invalid_return;

    
    auto get_components = [&](double tx, double ty, double tz) -> std::pair<double, double> {
    double phi0 = std::atan2(ty - circle_center.Y(), tx - circle_center.X());

    // Match the "winding" of the helix to the target roughly
    while (phi0 - last_angle >  TMath::Pi()) phi0 -= 2*TMath::Pi();
    while (phi0 - last_angle < -TMath::Pi()) phi0 += 2*TMath::Pi();

    // Helper lambda to calculate the exact 3D distance squared
    auto get_dist2 = [&](double p) {
        double hX = circle_center.X() + R_pred_cm * std::cos(p);
        double hY = circle_center.Y() + R_pred_cm * std::sin(p);
        double hZ = a * p + b;
        return std::pow(tx - hX, 2) + std::pow(ty - hY, 2) + std::pow(tz - hZ, 2);
    };

    // 1. GLOBAL CHECK: Ensure we are starting on the absolute closest coil
    // (Fixes cases where the target is physically closer to the coil above/below)
    double p_fine = phi0;
    double current_dist2 = get_dist2(p_fine);
    
 /*   for (double offset : {-2 * TMath::Pi(), 2 * TMath::Pi()}) {
        double test_p = phi0 + offset;
        double test_dist2 = get_dist2(test_p);
        if (test_dist2 < current_dist2) {
            current_dist2 = test_dist2;
            p_fine = test_p;
        }
    }*/

    // 2. LOCAL REFINEMENT: Dynamic 'while' loop instead of fixed 'for' loop
    int iter = 0;
    int max_iters = 15; // Safety limit to prevent infinite loops

    while (iter < max_iters) {
        double c = std::cos(p_fine), s = std::sin(p_fine);
        
        // Derivative of 3D distance squared
        double f = (tx - circle_center.X() - R_pred_cm*c)*(R_pred_cm*s) +
                   (ty - circle_center.Y() - R_pred_cm*s)*(-R_pred_cm*c) +
                   a*(a*p_fine + b - tz);
        double fp = R_pred_cm*R_pred_cm + a*a;
        
        // Try the new step
        double next_p = p_fine - (f / fp);
        double next_dist2 = get_dist2(next_p);

        // STOPPING CONDITION: 
        // If the step makes the distance bigger, or if the step is microscopic, we hit the minimum!
        if (next_dist2 >= current_dist2 || std::abs(next_p - p_fine) < 1e-6) {
            break; 
        }

        // Otherwise, accept the lower distance and continue
        p_fine = next_p;
        current_dist2 = next_dist2;
        iter++;
    }

    // Helix coordinates at the solved, verified minimum phi
        double hX = circle_center.X() + R_pred_cm * std::cos(p_fine);
        double hY = circle_center.Y() + R_pred_cm * std::sin(p_fine);
        double hZ = a * p_fine + b;

    // Calculate components separately
        double dXY = std::sqrt(std::pow(tx - hX, 2) + std::pow(ty - hY, 2));
        double dZ  = std::abs(tz - hZ);
    
        return {dXY, dZ};
    };
    
    
    // Helper to project to helix and get separated XY and Z components
   /* auto get_components = [&](double tx, double ty, double tz) -> std::pair<double, double> {
        double phi0 = std::atan2(ty - circle_center.Y(), tx - circle_center.X());
        
        // Match the "winding" of the helix to the target
        while (phi0 - last_angle >  TMath::Pi()) phi0 -= 2*TMath::Pi();
        while (phi0 - last_angle < -TMath::Pi()) phi0 += 2*TMath::Pi();

        // Newton solver to find the exact closest phi on the helix
        double p_fine = phi0;
        for (int iter = 0; iter < 5; ++iter) {
            double c = std::cos(p_fine), s = std::sin(p_fine);
            // Derivative of 3D distance squared
            double f = (tx - circle_center.X() - R_pred_cm*c)*(R_pred_cm*s) + 
                       (ty - circle_center.Y() - R_pred_cm*s)*(-R_pred_cm*c) + 
                       a*(a*p_fine + b - tz);
            double fp = R_pred_cm*R_pred_cm + a*a;
            p_fine -= f / fp;
        }


        // Helix coordinates at the solved phi
        double hX = circle_center.X() + R_pred_cm * std::cos(p_fine);
        double hY = circle_center.Y() + R_pred_cm * std::sin(p_fine);
        double hZ = a * p_fine + b;

        // Calculate components separately
        double dXY = std::sqrt(std::pow(tx - hX, 2) + std::pow(ty - hY, 2));
        double dZ  = std::abs(tz - hZ); 

        return {dXY, dZ};
    };*/

    // Calculate for MD0
    auto comp0 = get_components(data.md_anchor_x->at(md0), 
                                data.md_anchor_y->at(md0), 
                                data.md_anchor_z->at(md0));

    // Calculate for MD1
    auto comp1 = get_components(data.md_anchor_x->at(md1), 
                                data.md_anchor_y->at(md1), 
                                data.md_anchor_z->at(md1));

    // Return the 4-component vector
    return {comp0.first, comp0.second, comp1.first, comp1.second};
}

    std::pair<double, double> extrapolatePlsInRZAndGetDeltaR(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -999.0;
        const std::pair<double, double> invalid_pair = {-999.0, -999.0};

        std::vector<double> pls_z, pls_r;
        auto add_pls_hit = [&](float x, float y, float z) {
            if (x > invalid_val) {
                pls_z.push_back(z);
                pls_r.push_back(std::sqrt(x*x + y*y));
            }
        };
        add_pls_hit(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx));
        add_pls_hit(data.pls_hit1_x->at(pls_idx), data.pls_hit1_y->at(pls_idx), data.pls_hit1_z->at(pls_idx));
        add_pls_hit(data.pls_hit2_x->at(pls_idx), data.pls_hit2_y->at(pls_idx), data.pls_hit2_z->at(pls_idx));
        add_pls_hit(data.pls_hit3_x->at(pls_idx), data.pls_hit3_y->at(pls_idx), data.pls_hit3_z->at(pls_idx));

        if (pls_z.size() < 2) return invalid_pair;

        double sumZ = 0, sumR = 0, sumZR = 0, sumZ2 = 0;
        for (size_t i = 0; i < pls_z.size(); ++i) {
            sumZ += pls_z[i];
            sumR += pls_r[i];
            sumZR += pls_z[i] * pls_r[i];
            sumZ2 += pls_z[i] * pls_z[i];
        }
        double n = pls_z.size();
        double denom = (n * sumZ2 - sumZ * sumZ);
        if (std::abs(denom) < 1e-12) return invalid_pair;

        double slope = (n * sumZR - sumZ * sumR) / denom;
        double intercept = (sumR - slope * sumZ) / n;

        int md0 = data.ls_mdIdx0->at(ls_idx);
        int md1 = data.ls_mdIdx1->at(ls_idx);
        if (md0 < 0 || md1 < 0) return invalid_pair;

        double z_ls0 = data.md_anchor_z->at(md0);
        double r_ls0_act = std::sqrt(std::pow(data.md_anchor_x->at(md0), 2) + std::pow(data.md_anchor_y->at(md0), 2));
        double z_ls1 = data.md_anchor_z->at(md1);
        double r_ls1_act = std::sqrt(std::pow(data.md_anchor_x->at(md1), 2) + std::pow(data.md_anchor_y->at(md1), 2));

        return {r_ls0_act - (slope * z_ls0 + intercept), r_ls1_act - (slope * z_ls1 + intercept)};
    }

    std::pair<double, double> extrapolateSimplePointingInRZ(int pls_idx, int ls_idx, const rootReader& data) {
        const float invalid_val = -999.0;
        const std::pair<double, double> invalid_pair = {-999.0, -999.0};

        double z_first = 0, r_first = 0, z_last = 0, r_last = 0;
        bool found_first = false, found_last = false;

        if (data.pls_hit0_x->at(pls_idx) > invalid_val) {
            z_first = data.pls_hit0_z->at(pls_idx);
            r_first = std::sqrt(std::pow(data.pls_hit0_x->at(pls_idx), 2) + std::pow(data.pls_hit0_y->at(pls_idx), 2));
            found_first = true;
        } else if (data.pls_hit1_x->at(pls_idx) > invalid_val) {
            z_first = data.pls_hit1_z->at(pls_idx);
            r_first = std::sqrt(std::pow(data.pls_hit1_x->at(pls_idx), 2) + std::pow(data.pls_hit1_y->at(pls_idx), 2));
            found_first = true;
        }

        if (data.pls_hit3_x->at(pls_idx) > invalid_val) {
            z_last = data.pls_hit3_z->at(pls_idx);
            r_last = std::sqrt(std::pow(data.pls_hit3_x->at(pls_idx), 2) + std::pow(data.pls_hit3_y->at(pls_idx), 2));
            found_last = true;
        } else if (data.pls_hit2_x->at(pls_idx) > invalid_val) {
            z_last = data.pls_hit2_z->at(pls_idx);
            r_last = std::sqrt(std::pow(data.pls_hit2_x->at(pls_idx), 2) + std::pow(data.pls_hit2_y->at(pls_idx), 2));
            found_last = true;
        }

        if (!found_first || !found_last || std::abs(z_last - z_first) < 1e-4) return invalid_pair;

        double m = (r_last - r_first) / (z_last - z_first);
        double b = r_first - m * z_first;

        int md0 = data.ls_mdIdx0->at(ls_idx), md1 = data.ls_mdIdx1->at(ls_idx);
        if (md0 < 0 || md1 < 0) return invalid_pair;

        double z0 = data.md_anchor_z->at(md0), r0 = std::sqrt(std::pow(data.md_anchor_x->at(md0), 2) + std::pow(data.md_anchor_y->at(md0), 2));
        double z1 = data.md_anchor_z->at(md1), r1 = std::sqrt(std::pow(data.md_anchor_x->at(md1), 2) + std::pow(data.md_anchor_y->at(md1), 2));

        return {r0 - (m * z0 + b), r1 - (m * z1 + b)};
    }

    double calculateDeltaAngle(int pls_idx, int ls_idx, const rootReader& data) {
        const float inv = -999.0;
        TVector3 vPixel;
        if (data.pls_hit0_x->at(pls_idx) > inv) {
            vPixel.SetXYZ(data.pls_hit0_x->at(pls_idx), data.pls_hit0_y->at(pls_idx), data.pls_hit0_z->at(pls_idx));
        } else {
            return -999.0;
        }
        int md1_idx = data.ls_mdIdx1->at(ls_idx);
        if (md1_idx < 0) return -999.0;
        TVector3 vStrip(data.md_anchor_x->at(md1_idx), data.md_anchor_y->at(md1_idx), data.md_anchor_z->at(md1_idx));
        return vPixel.Angle(vStrip); 
    }

} // End namespace extrapolation
