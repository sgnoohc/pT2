#ifndef EXTRA_CUTS_H
#define EXTRA_CUTS_H

#include "rootReader.h"

namespace extra_cuts {

    // Helper math functions to exactly replicate cms::alpakatools math
    double reducePhiRange(double x);
    double deltaPhiLST(double x1, double y1, double x2, double y2);

    // Core LST Kinematic Calculations
    double calculateLSTDPhi(int pls_idx, int ls_idx, const rootReader& data);

    std::vector<double> calculateLSTdBeta(int pls_idx, int ls_idx, const rootReader& data);

    double calculateLSTKinematicZResidual(int pls_idx, int ls_idx, const rootReader& data);
    double calculateLSTOriginZResidual(int pls_idx, int ls_idx, const rootReader& data);
    int getCategoryFromDetId(uint32_t detId);
    int getConnectionIndex(int c0, int c1);

    }

#endif // LST_MATH_H
