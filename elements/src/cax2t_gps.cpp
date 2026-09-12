#include "cax2t_gps.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gauss = 0.577350269189625764509148780501957456;

void validate_geometry(const Cax2tGpsGeometry& geometry) {
    const double thickness = geometry.radii[1] - geometry.radii[0];
    const double height = geometry.z_upper - geometry.z_lower;
    if (!std::isfinite(geometry.radii[0]) || !std::isfinite(geometry.radii[1]) || !std::isfinite(geometry.z_lower)
        || !std::isfinite(geometry.z_upper) || geometry.radii[0] < 0.0 || !(thickness > 0.0)
        || !std::isfinite(thickness) || !(height > 0.0) || !std::isfinite(height))
        throw std::invalid_argument("CAX2T_GPS requires increasing nonnegative radii and increasing finite axial ends");
    const double volume = pi * thickness * (geometry.radii[0] + geometry.radii[1]) * height;
    if (!std::isfinite(volume) || !(volume > 0.0))
        throw std::invalid_argument("CAX2T_GPS reference volume must be positive and finite");
}

void validate_deformed_geometry(const Cax2tGpsGeometry& geometry, const Cax2tGpsLocalValues& state) {
    const double inner = geometry.radii[0] + state[2], outer = geometry.radii[1] + state[3];
    const double height = geometry.z_upper - geometry.z_lower + state[5] - state[4];
    // The axis condition and endpoint orientation also apply to small strain,
    // matching the production state validator. Positive Gauss-point radii alone
    // do not exclude an inner endpoint that has crossed the symmetry axis.
    if ((geometry.radii[0] == 0.0 && state[2] != 0.0) || !std::isfinite(inner) || !std::isfinite(outer)
        || !std::isfinite(height) || inner < 0.0 || !(outer > inner) || !(height > 0.0))
        throw std::domain_error(
            "CAX2T_GPS deformed endpoints must preserve the axis, ordered radii and positive axial height");
}

void validate_input(const Cax2tGpsInput& input) {
    validate_geometry(input.geometry);
    if (input.strain_formulation != StrainFormulation::small && input.strain_formulation != StrainFormulation::finite)
        throw std::invalid_argument("CAX2T_GPS has an unsupported strain formulation");
    if (!std::isfinite(input.time) || !std::isfinite(input.time_step) || input.time_step < 0.0
        || !std::isfinite(input.volumetric_heat_source))
        throw std::invalid_argument("CAX2T_GPS time, nonnegative time step and heat source must be finite");
    if (input.committed_history && !(input.time_step > 0.0))
        throw std::invalid_argument("CAX2T_GPS material history requires a positive time step");
    if (input.include_thermal_time_term && !input.committed_history)
        throw std::invalid_argument("CAX2T_GPS heat capacity requires committed material history");
    for (std::size_t i = 0; i < input.state.size(); ++i) {
        if (!std::isfinite(input.state[i]) || (i < 2 && !(input.state[i] > 0.0)))
            throw std::domain_error("CAX2T_GPS trial values must be finite and temperatures positive");
        if (input.committed_history
            && (!std::isfinite(input.committed_state[i]) || (i < 2 && !(input.committed_state[i] > 0.0))))
            throw std::invalid_argument("CAX2T_GPS committed values must be finite and temperatures positive");
        if (input.strain_formulation == StrainFormulation::finite && i >= 2 && !std::isfinite(input.committed_state[i]))
            throw std::invalid_argument("CAX2T_GPS committed displacements must be finite");
    }
    validate_deformed_geometry(input.geometry, input.state);
    if (input.committed_history || input.strain_formulation == StrainFormulation::finite)
        validate_deformed_geometry(input.geometry, input.committed_state);
}
} // namespace

Cax2tGpsGeometry make_cax2t_gps_geometry(const std::array<double, 2>& radii, double z_lower, double z_upper) {
    const Cax2tGpsGeometry geometry{radii, z_lower, z_upper};
    validate_geometry(geometry);
    return geometry;
}

Cax2tGpsResult evaluate_cax2t_gps(const Cax2tGpsInput& input, ElementRequest request) {
    validate_input(input);
    Cax2tGpsResult result;
    const auto& geometry = input.geometry;
    const auto& state = input.state;
    const bool finite = input.strain_formulation == StrainFormulation::finite;
    const double thickness = geometry.radii[1] - geometry.radii[0];
    const double height = geometry.z_upper - geometry.z_lower;
    const std::array<double, 2> gradient = {-1.0 / thickness, 1.0 / thickness};
    const double temperature_gradient = gradient[0] * state[0] + gradient[1] * state[1];
    const double radial_strain = gradient[0] * state[2] + gradient[1] * state[3];
    const double axial_strain = (state[5] - state[4]) / height;
    const double axial_coordinate = geometry.z_lower + 0.5 * height;
    const double expansion_temperature = 0.5 * state[0] + 0.5 * state[1];
    const std::array<double, 2> stations = {-gauss, gauss};
    for (std::size_t q = 0; q < stations.size(); ++q) {
        const std::array<double, 2> shape = {0.5 * (1.0 - stations[q]), 0.5 * (1.0 + stations[q])};
        const double radius = shape[0] * geometry.radii[0] + shape[1] * geometry.radii[1];
        double measure = pi * radius * thickness * height;
        const double temperature = shape[0] * state[0] + shape[1] * state[1];
        const double radial_displacement = shape[0] * state[2] + shape[1] * state[3];
        std::array<double, 4> strain = {radial_strain, axial_strain, radial_displacement / radius, 0.0};
        const MaterialFunctionContext context{input.time, radius, 0.0, axial_coordinate};
        // All geometry chains are explicit functions of three point-local
        // stretches; only the existing width-five material tangent uses AD.
        const std::array<Cax2tGpsLocalValues, 4> reference_chain = {{
            {{0.0, 0.0, gradient[0], gradient[1], 0.0, 0.0}},
            {{0.0, 0.0, 0.0, 0.0, -1.0 / height, 1.0 / height}},
            {{0.0, 0.0, shape[0] / radius, shape[1] / radius, 0.0, 0.0}},
            {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        }};
        auto strain_chain = reference_chain;
        auto virtual_strain = reference_chain;
        std::array<double, 4> stretch{1.0, 1.0, 1.0, 1.0};
        Cax2tGpsLocalValues measure_derivative{};
        if (finite) {
            const auto& old = input.committed_state;
            const std::array<double, 3> committed_stretch{1.0 + gradient[0] * old[2] + gradient[1] * old[3],
                1.0 + (old[5] - old[4]) / height,
                1.0 + (shape[0] * old[2] + shape[1] * old[3]) / radius};
            double committed_measure = measure, midpoint_measure = measure;
            for (std::size_t component = 0; component < 3; ++component) {
                stretch[component] = 1.0 + strain[component];
                const double midpoint = 0.5 * (stretch[component] + committed_stretch[component]);
                if (!(stretch[component] > 0.0) || !std::isfinite(stretch[component])
                    || !(committed_stretch[component] > 0.0) || !std::isfinite(committed_stretch[component])
                    || !(midpoint > 0.0) || !std::isfinite(midpoint))
                    throw std::domain_error("CAX2T_GPS current, committed and midpoint stretches must remain positive");
                committed_measure *= committed_stretch[component];
                midpoint_measure *= midpoint;
                // With diagonal F there is no spin. This is the diagonal
                // Hughes-Winget midpoint increment, not an accumulated log strain.
                strain[component] = (stretch[component] - committed_stretch[component]) / midpoint;
                const double increment_derivative = committed_stretch[component] / (midpoint * midpoint);
                for (std::size_t column = 0; column < state.size(); ++column) {
                    strain_chain[component][column] *= increment_derivative;
                    virtual_strain[component][column] /= stretch[component];
                }
            }
            measure *= stretch[0] * stretch[1] * stretch[2];
            if (!std::isfinite(measure) || !(measure > 0.0) || !std::isfinite(committed_measure)
                || !(committed_measure > 0.0) || !std::isfinite(midpoint_measure) || !(midpoint_measure > 0.0)
                || !std::isfinite(radius * stretch[2]) || !std::isfinite(radius * committed_stretch[2])
                || !std::isfinite(radius * (0.5 * stretch[2] + 0.5 * committed_stretch[2])))
                throw std::domain_error(
                    "CAX2T_GPS current, committed and midpoint radii and volumes must be finite and positive");
            if (request.jacobian)
                for (std::size_t column = 0; column < state.size(); ++column)
                    for (std::size_t component = 0; component < 3; ++component)
                        measure_derivative[column] += measure * reference_chain[component][column] / stretch[component];
            if (input.committed_history) {
                const double old_expansion_temperature = 0.5 * old[0] + 0.5 * old[1];
                auto old_context = context;
                old_context.time -= input.time_step;
                const auto old_eigenstrain = input.material.eigenstrain_rz(old_expansion_temperature, old_context);
                const std::array<double, 4> imposed = {old_eigenstrain.rr.value(),
                    old_eigenstrain.zz.value(),
                    old_eigenstrain.hoop.value(),
                    old_eigenstrain.rz.value()};
                const auto& old_history = (*input.committed_history)[q];
                for (std::size_t component = 0; component < 4; ++component)
                    strain[component] += old_history.elastic_strain[component] + old_history.plastic_strain[component]
                                         + old_history.creep_strain[component] + imposed[component];
            }
        }
        // The material interface retains point temperature for elastic and inelastic
        // properties. Correct its total-strain input so that it subtracts the
        // eigenstrain at the arithmetic mean of the two radial-node temperatures.
        // Only temperature is seeded here; the four-strain material seed stays width five.
        const adlite::Scalar point_temperature =
            request.jacobian ? adlite::Scalar::independent(temperature, 0, 1) : adlite::Scalar(temperature);
        const adlite::Scalar mean_temperature = request.jacobian
                                                    ? adlite::Scalar::independent(expansion_temperature, 0, 1)
                                                    : adlite::Scalar(expansion_temperature);
        const auto point_eigen = input.material.eigenstrain_rz(point_temperature, context);
        const auto mean_eigen = input.material.eigenstrain_rz(mean_temperature, context);
        const std::array<adlite::Scalar, 4> point_imposed = {point_eigen.rr,
            point_eigen.zz,
            point_eigen.hoop,
            point_eigen.rz};
        const std::array<adlite::Scalar, 4> mean_imposed = {mean_eigen.rr,
            mean_eigen.zz,
            mean_eigen.hoop,
            mean_eigen.rz};
        for (std::size_t component = 0; component < 4; ++component) {
            strain[component] += point_imposed[component].value() - mean_imposed[component].value();
            if (request.jacobian) {
                const double point_derivative =
                    point_imposed[component].derivative_size() == 0 ? 0.0 : point_imposed[component].derivative(0);
                const double mean_derivative =
                    mean_imposed[component].derivative_size() == 0 ? 0.0 : mean_imposed[component].derivative(0);
                for (std::size_t column = 0; column < 2; ++column)
                    strain_chain[component][column] += point_derivative * shape[column] - 0.5 * mean_derivative;
            }
        }
        const auto response = evaluate_axisymmetric_material_response(input.material,
            strain,
            temperature,
            input.time_step,
            input.committed_history ? &(*input.committed_history)[q] : nullptr,
            context,
            request.jacobian,
            request.history);
        if (request.history) {
            result.history[q] = response.history;
            result.history[q].stress = response.stress;
            if (!input.committed_history) {
                for (std::size_t component = 0; component < 4; ++component)
                    result.history[q].elastic_strain[component] = strain[component] - point_imposed[component].value();
            }
        }
        if (request.stress)
            result.stress[q] = response.stress;

        // The constant axial strain uses shared end-section values. It is
        // independent of radius, and introduces no radial-axial shear strain.
        const std::array<double, 4> stress = {response.stress.rr,
            response.stress.zz,
            response.stress.hoop,
            response.stress.rz};
        for (std::size_t row = 2; row < state.size(); ++row) {
            double force = 0.0;
            for (std::size_t component = 0; component < 4; ++component)
                force += virtual_strain[component][row] * stress[component];
            result.residual[row] += measure * force;
            if (request.jacobian)
                for (std::size_t column = 0; column < state.size(); ++column) {
                    double tangent = 0.0;
                    for (std::size_t component = 0; component < 4; ++component) {
                        double derivative = column < 2 ? response.thermal[component] * shape[column] : 0.0;
                        for (std::size_t strain_component = 0; strain_component < 4; ++strain_component)
                            derivative +=
                                response.tangent[component][strain_component] * strain_chain[strain_component][column];
                        tangent += virtual_strain[component][row] * derivative;
                        if (finite)
                            tangent -= virtual_strain[component][row] * reference_chain[component][column]
                                       / stretch[component] * stress[component];
                    }
                    result.jacobian[row * state.size() + column] +=
                        measure * tangent + measure_derivative[column] * force;
                }
        }

        // Only the integration-point temperature is seeded for thermal
        // properties; the two nodal derivatives are attached in closed form.
        const adlite::Scalar active_temperature =
            request.jacobian ? adlite::Scalar::independent(temperature, 0, 1) : adlite::Scalar(temperature);
        const adlite::Scalar conductivity = input.material.conductivity(active_temperature, context);
        const double conductivity_derivative =
            request.jacobian && conductivity.derivative_size() != 0 ? conductivity.derivative(0) : 0.0;
        double stored_rate = 0.0, stored_derivative = 0.0;
        if (input.include_thermal_time_term) {
            const double committed_temperature =
                shape[0] * input.committed_state[0] + shape[1] * input.committed_state[1];
            const double temperature_rate = (temperature - committed_temperature) / input.time_step;
            const adlite::Scalar capacity = input.material.heat_capacity(active_temperature, context);
            stored_rate = capacity.value() * temperature_rate;
            if (request.jacobian)
                stored_derivative =
                    capacity.value() / input.time_step
                    + (capacity.derivative_size() != 0 ? capacity.derivative(0) : 0.0) * temperature_rate;
        }
        result.stored_heat_rate += measure * stored_rate;
        result.generated_heat_rate += measure * input.volumetric_heat_source;
        const double current_temperature_gradient = temperature_gradient / stretch[0];
        for (std::size_t row = 0; row < 2; ++row) {
            const double current_gradient = gradient[row] / stretch[0];
            const double conduction = conductivity.value() * current_gradient * current_temperature_gradient;
            const double density = conduction + shape[row] * (stored_rate - input.volumetric_heat_source);
            result.residual[row] += measure * density;
            if (request.jacobian) {
                for (std::size_t column = 0; column < 2; ++column)
                    result.jacobian[row * state.size() + column] +=
                        measure
                        * (conductivity.value() * current_gradient * gradient[column] / stretch[0]
                            + conductivity_derivative * shape[column] * current_gradient * current_temperature_gradient
                            + shape[row] * stored_derivative * shape[column]);
                for (std::size_t column = 2; column < state.size(); ++column)
                    result.jacobian[row * state.size() + column] +=
                        measure_derivative[column] * density
                        - (finite ? 2.0 * measure * conduction * reference_chain[0][column] / stretch[0] : 0.0);
            }
        }
    }
    if (!request.residual && !request.jacobian)
        result.residual.fill(0.0);
    return result;
}
} // namespace fuelsim::elements
