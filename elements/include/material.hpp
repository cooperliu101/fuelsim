#pragma once
#include "element_types.hpp"
#include "material_types.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <memory>

namespace fuelsim {
struct ThermoelasticProperties final {
    std::shared_ptr<const MaterialFunctionSet> functions;
    double reference_young_modulus;
};

struct ActiveThermoelasticProperties final {
    adlite::Scalar lame_lambda, shear_modulus;
};

AxisymmetricStress rotate_axisymmetric_tensor(const AxisymmetricStress& tensor, const AxisymmetricRotation& rotation);
SymmetricTensor3 rotate_cartesian_tensor(const SymmetricTensor3& tensor, const CartesianRotation& rotation);

class IsotropicThermoelasticMaterial final {
  public:
    explicit IsotropicThermoelasticMaterial(ThermoelasticProperties properties);

    const MaterialFunctionSet& functions() const noexcept { return *_properties.functions; }

    adlite::Scalar conductivity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    adlite::Scalar heat_capacity(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    ActiveThermoelasticProperties active_properties(const adlite::Scalar& temperature,
        MaterialFunctionContext context = {}) const;
    AxisymmetricStrain eigenstrain_rz(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    SymmetricTensor3 eigenstrain(const adlite::Scalar& temperature, MaterialFunctionContext context = {}) const;
    AxisymmetricStress stress(const adlite::Scalar& strain_rr,
        const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop,
        const adlite::Scalar& strain_rz,
        const adlite::Scalar& temperature,
        MaterialFunctionContext context = {}) const;
    SymmetricTensor3 stress(const SymmetricTensor3& strain,
        const adlite::Scalar& temperature,
        MaterialFunctionContext context = {}) const;
    SymmetricTensor3Values
    stress_values(const SymmetricTensor3Values& strain, double temperature, MaterialFunctionContext context = {}) const;
    InelasticStressResponse response(const adlite::Scalar& strain_rr,
        const adlite::Scalar& strain_zz,
        const adlite::Scalar& strain_hoop,
        const adlite::Scalar& strain_rz,
        const adlite::Scalar& temperature,
        double time_step,
        const MaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianInelasticStressResponse response(const SymmetricTensor3& strain,
        const adlite::Scalar& temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianInelasticStressResponse incremental_response(const SymmetricTensor3& strain_increment,
        const CartesianRotation& rotation,
        const adlite::Scalar& temperature,
        double committed_temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianMaterialPointState response_values(const SymmetricTensor3Values& strain,
        double temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    CartesianMaterialPointState incremental_response_values(const SymmetricTensor3Values& strain_increment,
        const CartesianRotation& rotation,
        double temperature,
        double committed_temperature,
        double time_step,
        const CartesianMaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    InelasticStressResponse incremental_response(const adlite::Scalar& strain_increment_rr,
        const adlite::Scalar& strain_increment_zz,
        const adlite::Scalar& strain_increment_hoop,
        const adlite::Scalar& strain_increment_rz,
        const AxisymmetricRotation& rotation,
        const adlite::Scalar& temperature,
        double committed_temperature,
        double time_step,
        const MaterialPointState& committed,
        MaterialFunctionContext context = {}) const;
    static MaterialPointState state_values(const MaterialPointTrialState& trial_state);

  private:
    ThermoelasticProperties _properties;
};

struct AxisymmetricMaterialResponse final {
    AxisymmetricStressValues stress;
    std::array<std::array<double, 4>, 4> tangent{};
    std::array<double, 4> thermal{};
    // Unrotated trial strain histories; populated only when requested.
    MaterialPointState history{};
};

AxisymmetricMaterialResponse evaluate_axisymmetric_material_response(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 4>& fed_strain,
    double temperature,
    double time_step,
    const MaterialPointState* committed_material,
    MaterialFunctionContext context,
    bool compute_tangent,
    bool compute_history);

void rotate_axisymmetric_strain_history(MaterialPointState& history, const AxisymmetricRotation& rotation);

AxisymmetricStress compose_axisymmetric_stress(const AxisymmetricMaterialResponse& response,
    const std::array<adlite::Scalar, 5>& inputs,
    bool compute_tangent);

SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const std::array<std::array<double, 3>, 3>& rotation);
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
} // namespace fuelsim
