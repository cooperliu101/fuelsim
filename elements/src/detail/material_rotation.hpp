#pragma once
#include "cartesian_kinematics.hpp"
#include "material.hpp"

namespace fuelsim::element_detail {
SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const cartesian_detail::Matrix3& rotation);
CartesianInelasticStressResponse evaluate_incremental_cartesian_response(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain_increment,
    const adlite::Scalar& temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context);
} // namespace fuelsim::element_detail
