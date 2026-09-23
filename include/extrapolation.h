#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include <vector>
#include <utility>
#include <Math/Vector2D.h>
#include "rootReader.h" // Needed for the rootReader reference

namespace extrapolation {

    // Fits a circle of a known radius to a set of XY points to find the center (h, k)
    ROOT::Math::XYVector fitCircleWithFixedRadius(const std::vector<ROOT::Math::XYVector>& hits, double fixed_radius);

    // Main extrapolation function
    // Returns a pair representing the 3D distance to MD0 and MD1 respectively
    std::vector<double> extrapolatePlsHelicallyAndGetDistance(int pls_idx, int ls_idx, const rootReader& data);
    std::pair<double, double> extrapolatePlsInRZAndGetDeltaR(int pls_idx, int ls_idx, const rootReader& data);
    std::pair<double, double> extrapolateSimplePointingInRZ(int pls_idx, int ls_idx, const rootReader& data);
    double calculateDeltaAngle(int pls_idx, int ls_idx, const rootReader& data);

}

#endif // EXTRAPOLATION_H
