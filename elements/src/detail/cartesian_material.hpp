#pragma once
#include "coordinates.hpp"
#include "material.hpp"
#include "matrix3.hpp"

namespace fuelsim::element_detail {
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
} // namespace fuelsim::element_detail

namespace fuelsim::cartesian_detail {
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
