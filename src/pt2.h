#ifndef PT2_H
#define PT2_H

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

// A pLS + LS candidate and everything computed for it
struct pT2 {
    size_t pls_idx;
    size_t ls_idx;

    // Truth / kinematics
    float delta_pt = NAN;
    float delta_phi = NAN;
    float delta_eta = NAN;
    bool is_real = false;      // matched to sim?
    bool is_used = false;      // pLS or LS already used by an LST track candidate

    // Features, filled by computeFeatures()
    int combo_idx = -1;                     // layer-connection category, -1 if invalid
    int charge_idx = -1;                    // 0 = Pos, 1 = Neg
    float delta_r = NAN;
    float pls_eta = NAN, ls_eta = NAN;
    double delta_angle = NAN;
    std::vector<double> heli;               // MD0 dXY, MD0 dZ, MD1 dXY, MD1 dZ
    std::pair<double, double> rz_simple;    // MD0, MD1
    double lst_delta_phi = NAN, beta_in = NAN, beta_out = NAN, delta_beta = NAN, z_res_geo = NAN, z_res_kin = NAN;
    float nn_score = -1;                    // NN real-vs-fake score in [0, 1], -1 if no model loaded

    // Constructors
    pT2() = default;
    pT2(size_t pls, size_t ls) : pls_idx(pls), ls_idx(ls) {}

    void print(std::ostream &os = std::cout) const
    {
        os << "pT2 (pls " << pls_idx << ", ls " << ls_idx << "):\n"
            << "  delta_pt:      " << delta_pt << "\n"
            << "  delta_eta:     " << delta_eta << "\n"
            << "  delta_phi:     " << delta_phi << "\n"
            << "  is_real:       " << is_real << "\n"
            << "  is_used:       " << is_used << "\n"
            << "  combo_idx:     " << combo_idx << "\n"
            << "  charge_idx:    " << charge_idx << "\n"
            << "  delta_r:       " << delta_r << "\n"
            << "  pls_eta:       " << pls_eta << "\n"
            << "  ls_eta:        " << ls_eta << "\n"
            << "  delta_angle:   " << delta_angle << "\n"
            << "  heli:          ";
        for (double h : heli)
            os << " " << h;
        os << "\n"
            << "  rz_simple:     " << rz_simple.first << " " << rz_simple.second << "\n"
            << "  lst_delta_phi: " << lst_delta_phi << "\n"
            << "  beta_in:       " << beta_in << "\n"
            << "  delta_beta:    " << delta_beta << "\n"
            << "  beta_out:      " << beta_out << "\n"
            << "  z_res_geo:     " << z_res_geo << "\n"
            << "  z_res_kin:     " << z_res_kin << "\n"
            << "  nn_score:      " << nn_score << "\n";
    }
};

using pT2Collection = std::vector<pT2>;

#endif
