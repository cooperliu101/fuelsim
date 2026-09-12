#pragma once
#include "c3d20rt.hpp"
#include "c3d20t.hpp"
#include "c3d8rt.hpp"
#include "c3d8t.hpp"
#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "cax8rt.hpp"
#include "cax8t.hpp"
#include "contact_types.hpp"
#include "element_types.hpp"
#include "line2_rz.hpp"
#include "material.hpp"
#include "material_types.hpp"
#include <memory>
#include <stdexcept>
#include <utility>

namespace fuelsim::test {
inline bool same_material_state(const MaterialPointState& a, const MaterialPointState& b) {
    return a.elastic_strain == b.elastic_strain && a.plastic_strain == b.plastic_strain
           && a.creep_strain == b.creep_strain && a.equivalent_plastic_strain == b.equivalent_plastic_strain
           && a.equivalent_creep_strain == b.equivalent_creep_strain && a.stress.rr == b.stress.rr
           && a.stress.zz == b.stress.zz && a.stress.hoop == b.stress.hoop && a.stress.rz == b.stress.rz;
}

inline ThermoelasticProperties thermoelastic(double conductivity_inverse_temperature,
    double conductivity_offset,
    double young_modulus,
    double poisson_ratio,
    double thermal_expansion,
    double reference_temperature,
    double young_modulus_temperature_coefficient = 0.0,
    double poisson_ratio_temperature_coefficient = 0.0,
    double thermal_expansion_temperature_coefficient = 0.0,
    double density = 1.0,
    double specific_heat = 1.0) {
    MaterialFunctionRegistry registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = "test_material";
    functions->thermal = registry.bind_thermal("inverse_temperature_thermophysical",
        {{"conductivity_inverse_temperature", conductivity_inverse_temperature},
            {"conductivity_constant", conductivity_offset},
            {"density", density},
            {"specific_heat", specific_heat}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", young_modulus},
            {"poisson_ratio", poisson_ratio},
            {"reference_temperature", reference_temperature},
            {"young_modulus_temperature_coefficient", young_modulus_temperature_coefficient},
            {"poisson_ratio_temperature_coefficient", poisson_ratio_temperature_coefficient}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "linear_temperature_isotropic_thermal_expansion",
        {{"thermal_expansion", thermal_expansion},
            {"reference_temperature", reference_temperature},
            {"thermal_expansion_temperature_coefficient", thermal_expansion_temperature_coefficient}}));
    return {std::move(functions), young_modulus};
}

inline ThermoelasticProperties with_norton(ThermoelasticProperties material,
    double coefficient,
    double reference_stress,
    double stress_exponent,
    double reference_temperature = 600.0,
    double coefficient_temperature_coefficient = 0.0,
    double reference_stress_temperature_coefficient = 0.0,
    double stress_exponent_temperature_coefficient = 0.0) {
    MaterialFunctionRegistry registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>(*material.functions);
    functions->creep = registry.bind_creep("linear_temperature_norton",
        {{"coefficient", coefficient},
            {"reference_stress", reference_stress},
            {"stress_exponent", stress_exponent},
            {"reference_temperature", reference_temperature},
            {"coefficient_temperature_coefficient", coefficient_temperature_coefficient},
            {"reference_stress_temperature_coefficient", reference_stress_temperature_coefficient},
            {"stress_exponent_temperature_coefficient", stress_exponent_temperature_coefficient}});
    material.functions = std::move(functions);
    return material;
}

inline ThermoelasticProperties with_plasticity(ThermoelasticProperties material,
    double yield_stress,
    double hardening_modulus,
    double reference_temperature = 600.0,
    double yield_stress_temperature_coefficient = 0.0,
    double hardening_temperature_coefficient = 0.0) {
    MaterialFunctionRegistry registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>(*material.functions);
    functions->plasticity = registry.bind_plasticity("linear_temperature_isotropic_hardening",
        {{"yield_stress", yield_stress},
            {"hardening_modulus", hardening_modulus},
            {"reference_temperature", reference_temperature},
            {"yield_stress_temperature_coefficient", yield_stress_temperature_coefficient},
            {"hardening_temperature_coefficient", hardening_temperature_coefficient}});
    material.functions = std::move(functions);
    return material;
}
} // namespace fuelsim::test

namespace fuelsim {
// Test fixtures own materials and select the concrete model under test.
struct AxisymmetricTestData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source = 0.0, time = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    RzElementFormulation element_formulation = RzElementFormulation::cax4t;
    double initial_temperature = 600.0;
};

struct CartesianTestData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source, time;
    StrainFormulation strain_formulation = StrainFormulation::small;
    Hex8ElementFormulation hex8_element_formulation = Hex8ElementFormulation::c3d8t;
    double initial_temperature = 0.0;
};
} // namespace fuelsim

namespace fuelsim::rz {
struct LocalLinearization final {
    Cax4LocalResidual residual;
    Cax4LocalJacobian jacobian;
};

inline LocalLinearization linearize_line2_rz_boundary(const Line2RzBoundaryData& data,
    const Line2RzBoundaryGeometry& geometry,
    const Cax4LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_line2_rz_boundary(data, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_line2_rz_gap_heat(const GapHeatProperties& properties,
    const Line2RzHeatPointGeometry& geometry,
    const Cax4LocalValues& state) {
    LocalLinearization result{};
    result.residual = compute_line2_rz_gap_heat(properties, geometry, state, &result.jacobian);
    return result;
}

inline LocalLinearization linearize_node_to_line_rz_contact(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed_state,
    const ContactPointHistory& history) {
    LocalLinearization result{};
    result.residual =
        compute_node_to_line_rz_contact(properties, geometry, state, committed_state, history, &result.jacobian);
    return result;
}

} // namespace fuelsim::rz

namespace fuelsim {
inline elements::C3d8Result evaluate_c3d8(const CartesianTestData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    const CartesianMaterialHistory* history,
    double time_step,
    bool thermal_time,
    elements::ElementRequest request) {
    const elements::C3d8Input input{data.material,
        geometry,
        state,
        committed,
        history,
        time_step,
        data.time,
        data.volumetric_heat_source,
        data.strain_formulation,
        thermal_time,
        data.initial_temperature};
    return data.hex8_element_formulation == Hex8ElementFormulation::c3d8rt ? elements::evaluate_c3d8rt(input, request)
                                                                           : elements::evaluate_c3d8t(input, request);
}

inline Hex8LocalValues compute_c3d8_thermoelastic(const CartesianTestData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed = nullptr,
    double time_step = 0,
    Hex8LocalJacobian* jacobian = nullptr) {
    const auto result = evaluate_c3d8(data,
        geometry,
        state,
        committed ? *committed : Hex8LocalValues{},
        nullptr,
        time_step,
        time_step > 0,
        {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline Hex8LocalValues compute_c3d8_transient(const CartesianTestData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step,
    Hex8LocalJacobian* jacobian = nullptr,
    bool thermal_time = true) {
    const auto result = evaluate_c3d8(data,
        geometry,
        state,
        committed,
        &history,
        time_step,
        thermal_time,
        {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline CartesianMaterialHistory compute_c3d8_transient_update(const CartesianTestData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step) {
    return evaluate_c3d8(data, geometry, state, committed, &history, time_step, false, {false, false, true, false})
        .history;
}

inline std::array<SymmetricTensor3Values, 8>
compute_c3d8_stress(const CartesianTestData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state) {
    return evaluate_c3d8(data, geometry, state, {}, nullptr, 0, false, {false, false, false, true}).stress;
}

inline elements::C3d20Result evaluate_c3d20(const CartesianTestData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed,
    const CartesianMaterialHistory* history,
    double time_step,
    bool thermal_time,
    elements::ElementRequest request) {
    const elements::C3d20Input input{data.material,
        geometry,
        state,
        committed,
        history,
        time_step,
        data.time,
        data.volumetric_heat_source,
        data.strain_formulation,
        thermal_time,
        data.initial_temperature};
    return geometry.mechanical_points.size() == 8 ? elements::evaluate_c3d20rt(input, request)
                                                  : elements::evaluate_c3d20t(input, request);
}

inline Hex20LocalValues compute_c3d20_thermoelastic(const CartesianTestData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues* committed = nullptr,
    double time_step = 0,
    Hex20LocalJacobian* jacobian = nullptr) {
    const auto result = evaluate_c3d20(data,
        geometry,
        state,
        committed ? *committed : Hex20LocalValues{},
        nullptr,
        time_step,
        time_step > 0,
        {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline Hex20LocalValues compute_c3d20_transient(const CartesianTestData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step,
    Hex20LocalJacobian* jacobian = nullptr,
    bool thermal_time = true) {
    const auto result = evaluate_c3d20(data,
        geometry,
        state,
        committed,
        &history,
        time_step,
        thermal_time,
        {true, jacobian != nullptr, false, false});
    if (jacobian)
        *jacobian = result.jacobian;
    return result.residual;
}

inline CartesianMaterialHistory compute_c3d20_transient_update(const CartesianTestData& data,
    const Hex20Geometry& geometry,
    const Hex20LocalValues& state,
    const Hex20LocalValues& committed,
    const CartesianMaterialHistory& history,
    double time_step) {
    return evaluate_c3d20(data, geometry, state, committed, &history, time_step, false, {false, false, true, false})
        .history;
}

inline std::vector<SymmetricTensor3Values>
compute_c3d20_stress(const CartesianTestData& data, const Hex20Geometry& geometry, const Hex20LocalValues& state) {
    return evaluate_c3d20(data, geometry, state, {}, nullptr, 0, false, {false, false, false, true}).stress;
}
} // namespace fuelsim

namespace fuelsim::test {
inline Quad8RzGeometry make_cax8_geometry(const Quad8RzCoordinates& coordinates,
    RzElementFormulation model = RzElementFormulation::cax8t) {
    return model == RzElementFormulation::cax8rt ? elements::make_cax8rt_geometry(coordinates)
                                                 : elements::make_cax8t_geometry(coordinates);
}

inline Hex20Geometry make_c3d20_geometry(const Hex20Coordinates& coordinates,
    Hex20ElementFormulation model = Hex20ElementFormulation::c3d20t) {
    return model == Hex20ElementFormulation::c3d20rt ? elements::make_c3d20rt_geometry(coordinates)
                                                     : elements::make_c3d20t_geometry(coordinates);
}
} // namespace fuelsim::test
