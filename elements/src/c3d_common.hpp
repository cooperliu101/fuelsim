#pragma once
#include "c3d20_types.hpp"
#include "c3d8_types.hpp"
#include <adlite/adlite.hpp>
#include <array>

namespace fuelsim::cartesian_detail {
using Matrix3 = std::array<std::array<double, 3>, 3>;
using ActiveMatrix3 = std::array<std::array<adlite::Scalar, 3>, 3>;

double determinant(const Matrix3& matrix);
adlite::Scalar determinant(const ActiveMatrix3& matrix);
Matrix3 inverse(const Matrix3& matrix, double determinant_value);
ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& determinant_value);
ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second);

Matrix3 multiply(const Matrix3& first, const Matrix3& second);
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

namespace fuelsim::cartesian_detail {
SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const cartesian_detail::Matrix3& rotation);
SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const CartesianRotation& rotation);
CartesianInelasticStressResponse evaluate_incremental_cartesian_response(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain_increment,
    const adlite::Scalar& temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context);
MaterialFunctionContext material_context(double time, const CartesianPoint3& point);

struct CartesianStressTangent final {
    SymmetricTensor3Values stress;
    std::array<std::array<double, 6>, 6> tangent{};
    std::array<double, 6> thermal{};
};

CartesianStressTangent evaluate_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 6>& fed_strain,
    double temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context);
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

namespace fuelsim::c3d20_detail {
Hex20Geometry make_hex20_geometry(const Hex20Coordinates& coordinates, std::size_t order);
} // namespace fuelsim::c3d20_detail

namespace fuelsim::c3d20_detail {
elements::C3d20Result evaluate(const elements::C3d20Input& input, elements::ElementRequest request);
} // namespace fuelsim::c3d20_detail

namespace fuelsim::c3d8_detail {
void set_history_geometry(const elements::C3d8Input& input, bool reduced, elements::C3d8Result& result);
}

namespace fuelsim::c3d8_detail {
Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates);

void validate_cartesian_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
} // namespace fuelsim::c3d8_detail

namespace fuelsim::c3d20_detail {
void validate_hex20_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state);
} // namespace fuelsim::c3d20_detail
