#ifndef PT2_H
#define PT2_H

#include <vector>
#include <cmath>

struct pT2 {
    size_t pls_idx;
    size_t ls_idx;
    float delta_pt = NAN;
    float delta_phi = NAN;
    float delta_eta = NAN;    // eta difference
    bool is_real = false;      // matched to sim?
    bool is_used = false;      // for future processing

    // Constructors
    pT2() = default;
    pT2(size_t pls, size_t ls) : pls_idx(pls), ls_idx(ls) {}
};

using pT2Collection = std::vector<pT2>;

#endif

