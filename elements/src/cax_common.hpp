#pragma once
#include "cax4_types.hpp"

namespace fuelsim::cax4_detail {
Quad4RzGeometry make_quad4_rz_geometry(const Quad4Coordinates& coordinates);
} // namespace fuelsim::cax4_detail

namespace fuelsim {
struct AxisymmetricHughesWinget final {
    AxisymmetricRotation rotation;
    adlite::Scalar rr, zz, rz;
};

// Only the in-plane map; hoop strain and model-specific averaging remain with the caller.
AxisymmetricHughesWinget evaluate_axisymmetric_hughes_winget(const adlite::Scalar& hrr,
    const adlite::Scalar& hrz,
    const adlite::Scalar& hzr,
    const adlite::Scalar& hzz);
} // namespace fuelsim

namespace fuelsim {
struct AxisymmetricMidpointIncrement final {
    AxisymmetricHughesWinget in_plane;
    adlite::Scalar determinant_sum, hoop;
};

// Sum and difference are formed by the caller in its native variables, preserving
// the displacement-based or deformation-based arithmetic and averaging rules.
AxisymmetricMidpointIncrement evaluate_axisymmetric_midpoint_increment(const std::array<adlite::Scalar, 4>& sum,
    const std::array<adlite::Scalar, 4>& difference,
    const adlite::Scalar& hoop_sum,
    const adlite::Scalar& hoop_difference);
void finish_cax4_result(elements::Cax4Result& result, elements::ElementRequest request);
} // namespace fuelsim
