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
