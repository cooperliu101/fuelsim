#include "cartesian_material.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::cartesian_detail {
SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const cartesian_detail::Matrix3& rotation) {
    const cartesian_detail::Matrix3 value = {{{{tensor.xx, tensor.xy, tensor.xz}},
        {{tensor.xy, tensor.yy, tensor.yz}},
        {{tensor.xz, tensor.yz, tensor.zz}}}};
    const auto left = cartesian_detail::multiply(rotation, value);
    cartesian_detail::Matrix3 rotated{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                rotated[i][j] += left[i][k] * rotation[j][k];
    return {rotated[0][0], rotated[1][1], rotated[2][2], rotated[0][1], rotated[1][2], rotated[0][2]};
}

CartesianInelasticStressResponse evaluate_incremental_cartesian_response(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain_increment,
    const adlite::Scalar& temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) {
    if (!std::isfinite(committed_temperature) || !(committed_temperature > 0.0))
        throw std::domain_error("Incremental Cartesian material committed temperature must be finite and positive");
    MaterialFunctionContext old_context = context;
    old_context.time -= time_step;
    const SymmetricTensor3 old_imposed = material.eigenstrain(adlite::Scalar(committed_temperature), old_context);
    const SymmetricTensor3 synthetic_total{committed.elastic_strain[0] + strain_increment.xx + old_imposed.xx
                                               + committed.plastic_strain[0] + committed.creep_strain[0],
        committed.elastic_strain[1] + strain_increment.yy + old_imposed.yy + committed.plastic_strain[1]
            + committed.creep_strain[1],
        committed.elastic_strain[2] + strain_increment.zz + old_imposed.zz + committed.plastic_strain[2]
            + committed.creep_strain[2],
        committed.elastic_strain[3] + strain_increment.xy + old_imposed.xy + committed.plastic_strain[3]
            + committed.creep_strain[3],
        committed.elastic_strain[4] + strain_increment.yz + old_imposed.yz + committed.plastic_strain[4]
            + committed.creep_strain[4],
        committed.elastic_strain[5] + strain_increment.xz + old_imposed.xz + committed.plastic_strain[5]
            + committed.creep_strain[5]};
    return material.response(synthetic_total, temperature, time_step, committed, context);
}
} // namespace fuelsim::cartesian_detail

namespace fuelsim {
using namespace cartesian_detail;

SymmetricTensor3 rotate_cartesian_tensor(const SymmetricTensor3& tensor, const CartesianRotation& rotation) {
    const cartesian_detail::ActiveMatrix3 r = {{{rotation.xx, rotation.xy, rotation.xz},
        {rotation.yx, rotation.yy, rotation.yz},
        {rotation.zx, rotation.zy, rotation.zz}}};
    const cartesian_detail::ActiveMatrix3 value = {
        {{tensor.xx, tensor.xy, tensor.xz}, {tensor.xy, tensor.yy, tensor.yz}, {tensor.xz, tensor.yz, tensor.zz}}};
    cartesian_detail::ActiveMatrix3 left{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t l = 0; l < 3; ++l)
                left[i][k] += r[i][l] * value[l][k];
    cartesian_detail::ActiveMatrix3 rotated{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                rotated[i][j] += left[i][k] * r[j][k];
    return {rotated[0][0], rotated[1][1], rotated[2][2], rotated[0][1], rotated[1][2], rotated[0][2]};
}

CartesianInelasticStressResponse IsotropicThermoelasticMaterial::incremental_response(
    const SymmetricTensor3& strain_increment,
    const CartesianRotation& rotation,
    const adlite::Scalar& temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) const {
    CartesianInelasticStressResponse result = evaluate_incremental_cartesian_response(*this,
        strain_increment,
        temperature,
        committed_temperature,
        time_step,
        committed,
        context);
    result.stress = rotate_cartesian_tensor(result.stress, rotation);
    std::array<double, 6>* histories[3] = {&result.trial_state.elastic_strain,
        &result.trial_state.plastic_strain,
        &result.trial_state.creep_strain};
    for (std::array<double, 6>* history : histories) {
        const SymmetricTensor3 rotated = rotate_cartesian_tensor(
            {(*history)[0], (*history)[1], (*history)[2], (*history)[3], (*history)[4], (*history)[5]},
            rotation);
        *history = {rotated.xx.value(),
            rotated.yy.value(),
            rotated.zz.value(),
            rotated.xy.value(),
            rotated.yz.value(),
            rotated.xz.value()};
        for (double component : *history)
            if (!std::isfinite(component))
                throw std::domain_error("Rotated Cartesian material history must contain only finite values");
    }
    result.trial_state.stress = {result.stress.xx.value(),
        result.stress.yy.value(),
        result.stress.zz.value(),
        result.stress.xy.value(),
        result.stress.yz.value(),
        result.stress.xz.value()};
    return result;
}

CartesianMaterialPointState IsotropicThermoelasticMaterial::incremental_response_values(
    const SymmetricTensor3Values& strain_increment,
    const CartesianRotation& rotation,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) const {
    if (!std::isfinite(committed_temperature) || !(committed_temperature > 0.0))
        throw std::domain_error("Incremental Cartesian material committed temperature must be finite and positive");
    MaterialFunctionContext old_context = context;
    old_context.time -= time_step;
    const SymmetricTensor3 old_imposed = eigenstrain(adlite::Scalar(committed_temperature), old_context);
    const SymmetricTensor3Values synthetic_total{committed.elastic_strain[0] + strain_increment.xx
                                                     + old_imposed.xx.value() + committed.plastic_strain[0]
                                                     + committed.creep_strain[0],
        committed.elastic_strain[1] + strain_increment.yy + old_imposed.yy.value() + committed.plastic_strain[1]
            + committed.creep_strain[1],
        committed.elastic_strain[2] + strain_increment.zz + old_imposed.zz.value() + committed.plastic_strain[2]
            + committed.creep_strain[2],
        committed.elastic_strain[3] + strain_increment.xy + old_imposed.xy.value() + committed.plastic_strain[3]
            + committed.creep_strain[3],
        committed.elastic_strain[4] + strain_increment.yz + old_imposed.yz.value() + committed.plastic_strain[4]
            + committed.creep_strain[4],
        committed.elastic_strain[5] + strain_increment.xz + old_imposed.xz.value() + committed.plastic_strain[5]
            + committed.creep_strain[5]};
    CartesianMaterialPointState result = response_values(synthetic_total, temperature, time_step, committed, context);
    const cartesian_detail::Matrix3 rotation_values = {
        {{{rotation.xx.value(), rotation.xy.value(), rotation.xz.value()}},
            {{rotation.yx.value(), rotation.yy.value(), rotation.yz.value()}},
            {{rotation.zx.value(), rotation.zy.value(), rotation.zz.value()}}}};
    result.stress = rotate_cartesian_tensor_values(result.stress, rotation_values);
    std::array<double, 6>* histories[3] = {&result.elastic_strain, &result.plastic_strain, &result.creep_strain};
    for (std::array<double, 6>* history : histories) {
        const SymmetricTensor3Values rotated = rotate_cartesian_tensor_values(
            {(*history)[0], (*history)[1], (*history)[2], (*history)[3], (*history)[4], (*history)[5]},
            rotation_values);
        *history = {rotated.xx, rotated.yy, rotated.zz, rotated.xy, rotated.yz, rotated.xz};
        for (double component : *history)
            if (!std::isfinite(component))
                throw std::domain_error("Rotated Cartesian material history must contain only finite values");
    }
    return result;
}
} // namespace fuelsim

namespace fuelsim::cartesian_detail {
MaterialFunctionContext material_context(double time, const CartesianPoint3& point) {
    return {time, point.x, point.y, point.z};
}

CartesianStressTangent evaluate_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 6>& fed_strain,
    double temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context) {
    std::array<double, 7> seeds{};
    for (std::size_t component = 0; component < 6; ++component)
        seeds[component] = fed_strain[component];
    seeds[6] = temperature;
    std::array<adlite::Scalar, 7> active{};
    adlite::seed_identity(seeds.data(), seeds.size(), active.data());
    const SymmetricTensor3 strain{active[0], active[1], active[2], active[3], active[4], active[5]};
    const SymmetricTensor3 stress =
        committed_material == nullptr
            ? material.stress(strain, active[6], context)
            : material.response(strain, active[6], time_step, *committed_material, context).stress;
    const std::array<const adlite::Scalar*, 6> components =
        {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz};
    CartesianStressTangent result{};
    result.stress = {stress.xx.value(),
        stress.yy.value(),
        stress.zz.value(),
        stress.xy.value(),
        stress.yz.value(),
        stress.xz.value()};
    std::array<double, 7> derivatives{};
    for (std::size_t row = 0; row < 6; ++row) {
        components[row]->copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t column = 0; column < 6; ++column)
            result.tangent[row][column] = derivatives[column];
        result.thermal[row] = derivatives[6];
    }
    return result;
}

SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const CartesianRotation& rotation) {
    const cartesian_detail::Matrix3 values = {{{{rotation.xx.value(), rotation.xy.value(), rotation.xz.value()}},
        {{rotation.yx.value(), rotation.yy.value(), rotation.yz.value()}},
        {{rotation.zx.value(), rotation.zy.value(), rotation.zz.value()}}}};
    return rotate_cartesian_tensor_values(tensor, values);
}
} // namespace fuelsim::cartesian_detail
