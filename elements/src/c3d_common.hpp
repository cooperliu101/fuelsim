#pragma once
#include "c3d8_types.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace fuelsim::cartesian_detail {
using Matrix3 = std::array<std::array<double, 3>, 3>;
using ActiveMatrix3 = std::array<std::array<adlite::Scalar, 3>, 3>;

double determinant(const Matrix3& matrix);
adlite::Scalar determinant(const ActiveMatrix3& matrix);
Matrix3 inverse(const Matrix3& matrix, double determinant_value);
ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& determinant_value);
ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second);

Matrix3 multiply(const Matrix3& first, const Matrix3& second);

// The only project template exception: one Hughes-Winget map for primal and AD paths.
// Inputs are the central displacement gradient; geometry validation belongs to the caller.
template <typename Scalar>
void hughes_winget_rotation(const std::array<std::array<Scalar, 3>, 3>& gradient,
    std::array<std::array<Scalar, 3>, 3>& rotation,
    std::array<Scalar, 6>& strain) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>,
        "Hughes-Winget only supports double and adlite::Scalar");
    std::array<std::array<Scalar, 3>, 3> spatial_strain{}, numerator{}, denominator{};
    for (std::size_t i = 0; i < 3; ++i) {
        numerator[i][i] = 1.0;
        denominator[i][i] = 1.0;
        for (std::size_t j = 0; j < 3; ++j) {
            spatial_strain[i][j] = 0.5 * (gradient[i][j] + gradient[j][i]);
            const Scalar half_spin = 0.25 * (gradient[i][j] - gradient[j][i]);
            numerator[i][j] += half_spin;
            denominator[i][j] -= half_spin;
        }
    }
    const Scalar denominator_determinant = determinant(denominator);
    using std::isfinite;
    if (!isfinite(denominator_determinant) || denominator_determinant == 0.0)
        throw std::domain_error("Abaqus Hughes-Winget Cartesian rotation denominator is singular");
    const auto denominator_inverse = inverse(denominator, denominator_determinant);
    rotation = {};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                rotation[i][j] += numerator[i][k] * denominator_inverse[k][j];
    std::array<std::array<Scalar, 3>, 3> spatial_times_rotation{}, corotational_strain{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                spatial_times_rotation[i][j] += spatial_strain[i][k] * rotation[k][j];
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                corotational_strain[i][j] += rotation[k][i] * spatial_times_rotation[k][j];
    strain = {corotational_strain[0][0],
        corotational_strain[1][1],
        corotational_strain[2][2],
        corotational_strain[0][1],
        corotational_strain[1][2],
        corotational_strain[0][2]};
}
} // namespace fuelsim::cartesian_detail

namespace fuelsim::cartesian_detail {
struct KinematicsCore final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    ActiveMatrix3 current_inverse{};
    ActiveMatrix3 midpoint_inverse{};
    adlite::Scalar current_determinant{1.0};
};

KinematicsCore evaluate_kinematics(const ActiveMatrix3& gradient,
    const Matrix3& committed_deformation,
    StrainFormulation strain_formulation);
// This increment-only helper populates strain_increment and rotation. It cannot reconstruct the current
// configuration, so current_inverse and current_determinant retain their default values and must not be read.
KinematicsCore evaluate_hughes_winget_increment(const ActiveMatrix3& central_displacement_gradient);

} // namespace fuelsim::cartesian_detail

namespace fuelsim::c3d8_detail {
constexpr std::array<std::array<double, 3>, 8> hex8_signs = {{{{-1.0, -1.0, -1.0}},
    {{1.0, -1.0, -1.0}},
    {{1.0, 1.0, -1.0}},
    {{-1.0, 1.0, -1.0}},
    {{-1.0, -1.0, 1.0}},
    {{1.0, -1.0, 1.0}},
    {{1.0, 1.0, 1.0}},
    {{-1.0, 1.0, 1.0}}}};
// This permutation is its own inverse: node order to Gauss order and back.
constexpr std::array<std::size_t, hex8_node_count> hex8_node_gauss_permutation = {0, 1, 3, 2, 4, 5, 7, 6};

} // namespace fuelsim::c3d8_detail

namespace fuelsim::c3d8_detail {
struct C3d8Kinematics final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> current_gradient;
    adlite::Scalar current_weighted_measure{0.0};
};

C3d8Kinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation);

cartesian_detail::ActiveMatrix3 displacement_gradient(const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state);
cartesian_detail::Matrix3 deformation_gradient(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
C3d8Kinematics evaluate_cartesian_kinematics_from_gradient(const Hex8QuadraturePoint& point,
    const cartesian_detail::ActiveMatrix3& gradient,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation);
} // namespace fuelsim::c3d8_detail

namespace fuelsim::c3d8_detail {
void set_history_geometry(const elements::C3d8Input& input, bool reduced, elements::C3d8Result& result);
}

namespace fuelsim::c3d8_detail {
Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates);

void validate_cartesian_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
} // namespace fuelsim::c3d8_detail
