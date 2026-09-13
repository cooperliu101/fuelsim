#include "c3d8rt.hpp"
#include "support/c3d8_common_tests.hpp"

namespace {
using namespace fuelsim::test::c3d8;

bool test_reduced_integration_inelastic_jacobian() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8rt_geometry(coordinates);
    fuelsim::Hex8LocalValues committed_state{}, state{};
    for (std::size_t node = 0; node < 8; ++node) {
        committed_state[node] = 300.0;
        state[node] = 302.0 + 0.35 * static_cast<double>(node);
        const fuelsim::CartesianPoint3& point = coordinates[node];
        state[8 + node] = 0.20 * point.x + 0.02 * point.y + 0.01 * point.z;
        state[16 + node] = 0.02 * point.x - 0.04 * point.y - 0.015 * point.z;
        state[24 + node] = 0.01 * point.x - 0.015 * point.y - 0.03 * point.z;
    }
    const fuelsim::CartesianMaterialHistory committed_material(1);
    bool passed = true;
    for (const std::array<bool, 2> branch : {std::array<bool, 2>{false, true}, {true, false}, {true, true}}) {
        const fuelsim::CartesianTestData data{
            fuelsim::IsotropicThermoelasticMaterial(inelastic_properties(branch[0], branch[1])),
            3.0,
            1.0,
            fuelsim::StrainFormulation::small,
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        fuelsim::Hex8LocalJacobian jacobian{};
        (void)
            fuelsim::compute_c3d8_transient(data, geometry, state, committed_state, committed_material, 1.0, &jacobian);

        std::array<double, 32> direction{};
        for (std::size_t dof = 0; dof < direction.size(); ++dof)
            direction[dof] = std::sin(0.31 * static_cast<double>(dof + 1));
        std::array<double, 3> block_errors{};
        for (std::size_t block = 0; block < block_errors.size(); ++block) {
            const bool thermal_columns = block != 2;
            const bool thermal_rows = block == 0;
            const double step = thermal_columns ? 1.0e-4 : 1.0e-7;
            fuelsim::Hex8LocalValues plus = state, minus = state;
            for (std::size_t column = 0; column < direction.size(); ++column) {
                if ((column < 8) != thermal_columns)
                    continue;
                plus[column] += step * direction[column];
                minus[column] -= step * direction[column];
            }
            const fuelsim::Hex8LocalResidual plus_residual =
                fuelsim::compute_c3d8_transient(data, geometry, plus, committed_state, committed_material, 1.0);
            const fuelsim::Hex8LocalResidual minus_residual =
                fuelsim::compute_c3d8_transient(data, geometry, minus, committed_state, committed_material, 1.0);
            double difference_squared = 0.0, reference_squared = 0.0;
            for (std::size_t row = 0; row < 32; ++row) {
                if ((row < 8) != thermal_rows)
                    continue;
                double analytic = 0.0;
                for (std::size_t column = 0; column < 32; ++column) {
                    if ((column < 8) != thermal_columns)
                        continue;
                    analytic += jacobian[row * 32 + column] * direction[column];
                }
                const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
                const double difference = analytic - numerical;
                difference_squared += difference * difference;
                reference_squared += numerical * numerical;
            }
            block_errors[block] = std::sqrt(difference_squared / reference_squared);
        }
        double thermal_displacement_maximum = 0.0;
        for (std::size_t row = 0; row < 8; ++row)
            for (std::size_t column = 8; column < 32; ++column)
                thermal_displacement_maximum =
                    std::max(thermal_displacement_maximum, std::abs(jacobian[row * 32 + column]));
        const fuelsim::CartesianMaterialHistory update =
            fuelsim::compute_c3d8_transient_update(data, geometry, state, committed_state, committed_material, 1.0);
        const fuelsim::CartesianMaterialPointState& point = update.front();
        const double plastic_trace = point.plastic_strain[0] + point.plastic_strain[1] + point.plastic_strain[2];
        const double creep_trace = point.creep_strain[0] + point.creep_strain[1] + point.creep_strain[2];
        std::cout << "hex8_c3d8rt_inelastic_ktt_relative_error=" << block_errors[0] << '\n'
                  << "hex8_c3d8rt_inelastic_kut_relative_error=" << block_errors[1] << '\n'
                  << "hex8_c3d8rt_inelastic_kuu_relative_error=" << block_errors[2] << '\n';
        passed =
            check(*std::max_element(block_errors.begin(), block_errors.end()) < 2.0e-6,
                "C3D8RT plastic, creep, and coupled Jacobian blocks match centered directional differences")
            && check(thermal_displacement_maximum == 0.0,
                "small-strain C3D8RT thermal residual has no displacement coupling")
            && check(update.size() == 1
                         && (branch[1] ? point.equivalent_plastic_strain > 0.0 : point.equivalent_plastic_strain == 0.0)
                         && (branch[0] ? point.equivalent_creep_strain > 0.0 : point.equivalent_creep_strain == 0.0)
                         && std::abs(plastic_trace) < 2.0e-15 && std::abs(creep_trace) < 2.0e-15,
                "C3D8RT commits exactly one traceless plastic-creep material point with the requested branches")
            && check(committed_material.front().equivalent_plastic_strain == 0.0
                         && committed_material.front().equivalent_creep_strain == 0.0,
                "C3D8RT residual, Jacobian, and trial update do not mutate committed history")
            && passed;
    }
    const fuelsim::CartesianTestData finite_data{fuelsim::IsotropicThermoelasticMaterial(properties()),
        3.0,
        1.0,
        fuelsim::StrainFormulation::finite,
        fuelsim::Hex8ElementFormulation::c3d8rt,
        300.0};
    fuelsim::Hex8LocalJacobian finite_jacobian{};
    const fuelsim::Hex8LocalResidual finite_residual = fuelsim::compute_c3d8_transient(finite_data,
        geometry,
        state,
        committed_state,
        committed_material,
        1.0,
        &finite_jacobian);
    const fuelsim::Hex8LocalResidual finite_residual_without_jacobian =
        fuelsim::compute_c3d8_transient(finite_data, geometry, state, committed_state, committed_material, 1.0);
    const fuelsim::CartesianMaterialHistory finite_update =
        fuelsim::compute_c3d8_transient_update(finite_data, geometry, state, committed_state, committed_material, 1.0);
    double thermal_displacement_maximum = 0.0;
    for (std::size_t row = 0; row < 8; ++row)
        for (std::size_t column = 8; column < 32; ++column)
            thermal_displacement_maximum =
                std::max(thermal_displacement_maximum, std::abs(finite_jacobian[row * 32 + column]));
    passed = check(std::all_of(finite_residual.begin(),
                       finite_residual.end(),
                       [](double value) { return std::isfinite(value); })
                       && thermal_displacement_maximum > 0.0 && finite_update.size() == 1,
                 "finite-strain C3D8RT evaluates current-geometry thermal-mechanical coupling and one material point")
             && check(finite_residual == finite_residual_without_jacobian,
                 "finite-strain C3D8RT Jacobian and residual-only paths produce identical residual values")
             && passed;
    const auto rejects_domain = [&](const fuelsim::Hex8LocalValues& trial, const fuelsim::Hex8LocalValues& committed) {
        try {
            (void)fuelsim::compute_c3d8_transient(finite_data, geometry, trial, committed, committed_material, 1.0);
            return false;
        } catch (const std::domain_error&) {
            return true;
        }
    };
    fuelsim::Hex8LocalValues invalid_current = state;
    for (std::size_t node = 0; node < 8; ++node)
        invalid_current[8 + node] = -2.0 * coordinates[node].x;
    fuelsim::Hex8LocalValues invalid_midpoint = committed_state;
    for (std::size_t node = 0; node < 8; ++node) {
        invalid_midpoint[8 + node] = -2.0 * coordinates[node].x;
        invalid_midpoint[16 + node] = -2.0 * coordinates[node].y;
    }
    fuelsim::Hex8LocalValues invalid_committed = committed_state, expanded_current = committed_state;
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t component = 0; component < 3; ++component) {
            const double coordinate = component == 0   ? coordinates[node].x
                                      : component == 1 ? coordinates[node].y
                                                       : coordinates[node].z;
            invalid_committed[8 * (component + 1) + node] = -2.0 * coordinate;
            expanded_current[8 * (component + 1) + node] = 2.0 * coordinate;
        }
    passed = check(rejects_domain(invalid_current, committed_state) && rejects_domain(invalid_midpoint, committed_state)
                       && rejects_domain(expanded_current, invalid_committed),
                 "finite-strain C3D8RT rejects nonpositive current, midpoint, and committed configurations")
             && passed;
    return passed;
}

bool test_reduced_integration_thermoelastic_capacity_gate() {
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8rt_geometry(unit_cube());
    fuelsim::Hex8LocalValues committed_state{}, state{};
    for (std::size_t node = 0; node < 8; ++node) {
        committed_state[node] = 300.0;
        state[node] = 310.0;
    }
    std::array<fuelsim::Hex8LocalResidual, 2> steady{}, transient{};
    std::array<fuelsim::Hex8LocalJacobian, 2> jacobian{};
    const std::array<fuelsim::StrainFormulation, 2> formulations = {fuelsim::StrainFormulation::small,
        fuelsim::StrainFormulation::finite};
    bool passed = true;
    for (std::size_t formulation = 0; formulation < formulations.size(); ++formulation) {
        const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(capacity_properties()),
            0.0,
            1.0,
            formulations[formulation],
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        steady[formulation] = fuelsim::compute_c3d8_thermoelastic(data, geometry, state);
        const fuelsim::elements::C3d8Input input{data.material,
            geometry,
            state,
            committed_state,
            nullptr,
            2.0,
            data.time,
            0.0,
            formulations[formulation],
            false,
            300.0};
        const auto without_capacity = fuelsim::elements::evaluate_c3d8rt(input, {true, true, false, false});
        passed = check(without_capacity.residual == steady[formulation],
                     "C3D8RT model request can disable capacity while retaining the supplied committed configuration")
                 && passed;
        passed = check(without_capacity.history.empty(), "Unrequested C3D8RT history stays empty") && passed;

        transient[formulation] =
            fuelsim::compute_c3d8_thermoelastic(data, geometry, state, &committed_state, 2.0, &jacobian[formulation]);
        for (std::size_t node = 0; node < 8; ++node)
            passed = check(transient[formulation][node] - steady[formulation][node] > 0.0,
                         "C3D8RT thermoelastic local interface adds heat capacity when a committed state is supplied")
                     && passed;
    }
    for (std::size_t node = 0; node < 8; ++node) {
        passed = check(near(transient[0][node] - steady[0][node], transient[1][node] - steady[1][node], 1.0e-13),
                     "small- and finite-strain C3D8RT use the same committed-state heat-capacity gate")
                 && passed;
        for (std::size_t column = 0; column < 8; ++column)
            passed = check(near(jacobian[0][node * 32 + column], jacobian[1][node * 32 + column], 1.0e-13),
                         "undeformed small- and finite-strain C3D8RT thermal Jacobians include the same capacity term")
                     && passed;
    }
    return passed;
}

bool test_reduced_integration_hourglass_energy() {
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8rt_geometry(unit_cube());
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = 400.0;
        const double mode = (node == 0 || node == 2 || node == 5 || node == 7) ? 1.0 : -1.0;
        state[8 + node] = 0.03 * mode + 0.01 * geometry.reduced_point.gradient[node][0];
        state[16 + node] = -0.02 * mode;
        state[24 + node] = 0.015 * mode;
    }
    const fuelsim::ThermoelasticProperties fixed = fuelsim::test::thermoelastic(0.0, 1.0, 200.0, 0.25, 0.0, 300.0, 0.0);
    const fuelsim::ThermoelasticProperties varying_initial =
        fuelsim::test::thermoelastic(0.0, 1.0, 100.0, 0.25, 0.0, 300.0, 1.0);
    bool passed = true;
    for (const fuelsim::StrainFormulation formulation :
        {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        const fuelsim::CartesianTestData fixed_data{fuelsim::IsotropicThermoelasticMaterial(fixed),
            0.0,
            0.0,
            formulation,
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        const fuelsim::CartesianTestData varying_data{fuelsim::IsotropicThermoelasticMaterial(varying_initial),
            0.0,
            0.0,
            formulation,
            fuelsim::Hex8ElementFormulation::c3d8rt,
            300.0};
        const fuelsim::Hex8LocalResidual fixed_residual =
            fuelsim::compute_c3d8_thermoelastic(fixed_data, geometry, state);
        const fuelsim::Hex8LocalResidual varying_residual =
            fuelsim::compute_c3d8_thermoelastic(varying_data, geometry, state);
        double maximum_error = 0.0, maximum_scale = 0.0;
        constexpr double step = 1.0e-7;
        for (std::size_t column = 8; column < state.size(); ++column) {
            fuelsim::Hex8LocalValues plus = state, minus = state;
            plus[column] += step;
            minus[column] -= step;
            const double plus_energy = fuelsim::elements::c3d8rt_hourglass_energy({fixed_data.material,
                                           geometry,
                                           plus,
                                           plus,
                                           nullptr,
                                           0.0,
                                           fixed_data.time,
                                           fixed_data.volumetric_heat_source,
                                           fixed_data.strain_formulation,
                                           false,
                                           fixed_data.initial_temperature})
                                       - fuelsim::elements::c3d8rt_hourglass_energy({varying_data.material,
                                           geometry,
                                           plus,
                                           plus,
                                           nullptr,
                                           0.0,
                                           varying_data.time,
                                           varying_data.volumetric_heat_source,
                                           varying_data.strain_formulation,
                                           false,
                                           varying_data.initial_temperature});
            const double minus_energy = fuelsim::elements::c3d8rt_hourglass_energy({fixed_data.material,
                                            geometry,
                                            minus,
                                            minus,
                                            nullptr,
                                            0.0,
                                            fixed_data.time,
                                            fixed_data.volumetric_heat_source,
                                            fixed_data.strain_formulation,
                                            false,
                                            fixed_data.initial_temperature})
                                        - fuelsim::elements::c3d8rt_hourglass_energy({varying_data.material,
                                            geometry,
                                            minus,
                                            minus,
                                            nullptr,
                                            0.0,
                                            varying_data.time,
                                            varying_data.volumetric_heat_source,
                                            varying_data.strain_formulation,
                                            false,
                                            varying_data.initial_temperature});
            const double numerical = (plus_energy - minus_energy) / (2.0 * step);
            const double analytic = fixed_residual[column] - varying_residual[column];
            maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
            maximum_scale = std::max({maximum_scale, std::abs(analytic), std::abs(numerical)});
        }
        const double relative_error = maximum_error / maximum_scale;
        std::cout << "hex8_c3d8rt_hourglass_energy_gradient_relative_error=" << relative_error << '\n';
        passed = check(relative_error < 2.0e-8,
                     "C3D8RT mechanical hourglass energy gradient equals the residual hourglass force")
                 && passed;
    }
    return passed;
}

bool test_reduced_finite_native_increment(const fuelsim::Hex8Geometry& geometry) {
    // Element 1, increments 1 and 2 of the same native isothermal disconnected
    // fixture used below. The CSV exports S but not EE; both elastic histories
    // are independently obtained from S using constant isotropic compliance.
    const std::array<std::array<std::array<double, 3>, 8>, 2> displacement = {
        {{{{-4.3065502324946349e-32, -2.3315518925463993e-33, -4.8686329105615769e-33},
             {0.0049401908423474354, 0.0066122503072034685, -0.0034133679019782192},
             {-0.0099282070235988049, 0.006889888274352691, -0.003289889290425558},
             {-3.9171775421242442e-32, 2.3478032278779983e-33, -4.9880915392523464e-33},
             {-7.7481541307303104e-34, -1.7429875981271938e-33, -1.548059350628101e-33},
             {0.0061078998955587391, 0.0070002206920566573, -0.0026930759301446879},
             {-0.0061742457723674921, 0.0070674375322660474, -0.0033772986333970539},
             {-2.4485820501317906e-33, 8.4596564330026235e-34, -2.4200218183146214e-33}}},
            {{{3.7320354437479699e-32, -1.3425646017450348e-34, 9.0789508430568721e-35},
                {0.0049467520782551337, 0.0066086872361525544, -0.0034338577055195579},
                {-0.0099333491445455536, 0.0069032665498241367, -0.0033007786900270545},
                {3.4834029409267245e-32, 4.8889496068782318e-33, -2.1455222295399895e-33},
                {5.6611967043606205e-32, -2.5714639423299086e-33, 4.5254501209717765e-33},
                {0.0061208070715082566, 0.0069935763887165889, -0.0027026335427850399},
                {-0.0061653751946648188, 0.0070804406650137311, -0.0033798344638911182},
                {5.5232992305606835e-32, 4.0127638381778732e-34, 2.5021889049505611e-33}}}}};
    const std::array<std::array<double, 6>, 2> native_stress = {{{-281793.51171911135,
                                                                     -10230.861839389072,
                                                                     -4734.3699683770392,
                                                                     -2653.3650280368552,
                                                                     -250.54142425023213,
                                                                     -19257.332753339699},
        {-278910.73652076186,
            -6765.5215568881104,
            -1264.0646380595863,
            -2651.8545745333431,
            -174.01137881597074,
            -19284.704125947548}}};
    constexpr double young_modulus = 2.0e8, poisson_ratio = 0.25;
    std::array<fuelsim::Hex8LocalValues, 2> state{};
    std::array<std::array<double, 6>, 2> elastic_strain{};
    for (std::size_t frame = 0; frame < 2; ++frame) {
        const auto& stress = native_stress[frame];
        const double trace = stress[0] + stress[1] + stress[2];
        for (std::size_t component = 0; component < 6; ++component)
            elastic_strain[frame][component] =
                ((1.0 + poisson_ratio) * stress[component] - (component < 3 ? poisson_ratio * trace : 0.0))
                / young_modulus;
        for (std::size_t node = 0; node < 8; ++node) {
            state[frame][node] = 300.0;
            for (std::size_t component = 0; component < 3; ++component)
                state[frame][8 * (component + 1) + node] = displacement[frame][node][component];
        }
    }
    fuelsim::CartesianMaterialHistory history(1);
    history[0].elastic_strain = elastic_strain[0];
    const fuelsim::IsotropicThermoelasticMaterial material(
        fuelsim::test::thermoelastic(0.0, 1.0, young_modulus, poisson_ratio, 1.0e-5, 300.0));
    const fuelsim::elements::C3d8Input input{material,
        geometry,
        state[1],
        state[0],
        &history,
        100000.0,
        200000.0,
        0.0,
        fuelsim::StrainFormulation::finite,
        false,
        300.0};
    const auto result = fuelsim::elements::evaluate_c3d8rt(input, {false, false, true, false});
    if (!check(result.history.size() == 1, "native C3D8RT increment updates one active material point"))
        return false;
    const auto& updated = result.history[0];
    const auto& stress = updated.stress;
    const std::array<double, 6> actual_stress = {stress.xx, stress.yy, stress.zz, stress.xy, stress.yz, stress.xz};
    double strain_error = 0.0, stress_relative_error = 0.0;
    for (std::size_t component = 0; component < 6; ++component) {
        strain_error =
            std::max(strain_error, std::abs(updated.elastic_strain[component] - elastic_strain[1][component]));
        stress_relative_error = std::max(stress_relative_error,
            std::abs((actual_stress[component] - native_stress[1][component]) / native_stress[1][component]));
    }
    std::cout << "hex8_c3d8rt_native_increment_elastic_strain_absolute_error=" << strain_error << '\n'
              << "hex8_c3d8rt_native_increment_stress_maximum_relative_error=" << stress_relative_error << '\n';
    return check(strain_error < 1.0e-14 && stress_relative_error < 1.0e-9,
               "native nonaffine finite C3D8RT increment matches all six elastic strains and stresses")
           && check(history[0].elastic_strain == elastic_strain[0],
               "native C3D8RT increment leaves the supplied elastic history unchanged");
}

bool test_reduced_finite_native_bulk_force() {
    // Native fixture: verification/abaqus/b544_diagnosis/
    // b544_disconnected_weak_isothermal_{nodal,integration}.csv, increment 10, element 1.
    // All nodal temperatures are 300 K; E = 2e8 Pa, nu = .25. The eight nodes are
    // disconnected from adjacent elements, so their RF values are individual element forces.
    const fuelsim::Hex8Coordinates coordinates = {{{0.0, 0.0, 0.0},
        {1.06, -0.04, -0.025},
        {0.94, 0.96, 0.025},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 0.5},
        {1.0, 0.0, 0.475},
        {1.0, 1.0, 0.525},
        {0.0, 1.0, 0.5}}};
    const std::array<std::array<double, 3>, 8> displacement = {
        {{7.9462602585974087e-35, 8.6971022351024285e-35, 1.7069499058260809e-34},
            {0.0050666033739808734, 0.006414311154585366, -0.003751236305509646},
            {-0.0098916993558567975, 0.0071531224563867284, -0.0035800725755716972},
            {7.4166263286936459e-35, -8.4614140109446767e-35, 1.6761668268700088e-34},
            {3.1245610097098389e-32, 1.7483078954470366e-33, -1.866827126704909e-33},
            {0.0064339625085310876, 0.0068351920960517858, -0.0027764690216448618},
            {-0.0059642275370004744, 0.0073580788203811672, -0.0034634957136284489},
            {2.9242728979686527e-32, -6.6273233072827782e-34, -2.4519351027982266e-33}}};
    const std::array<double, 6> native_stress = {-202921.75001341777,
        77474.112214676948,
        82905.698002405756,
        -2774.7440110752673,
        405.38407577923766,
        -19908.488695006072};
    const std::array<std::array<double, 3>, 8> native_reaction = {
        {{31779.156494033879, -9550.0311541470455, -18229.148721612106},
            {-19127.317219542001, -9017.4576594611935, -22608.421222965255},
            {-21391.174521525558, 10314.606006157388, -23814.787661992275},
            {29172.080801707467, 9771.0188400530369, -18059.890982613339},
            {20249.321058793117, -9801.9778399960815, 23409.485626921964},
            {-30604.646347611189, -10301.25086330307, 19024.536510915277},
            {-30294.084004257886, 9048.4164099165773, 17258.826158089923},
            {20216.663738402167, 9536.6762607803976, 23019.400293255811}}};
    // The native card specifies an absolute weak hourglass stiffness of 1 Pa.
    // These corrections are fixed evaluations of the independently identified B533
    // hourglass energy derivative, via analyze_bulk_split.py::operators. Neither
    // these values nor the native RF fixture use a Fuelsim bulk-force calculation.
    const std::array<std::array<double, 3>, 8> weak_hourglass_force = {
        {{-0.001322365571272695, 0.00015753958571390346, 1.2611594021823463e-05},
            {0.0014893175092933813, -0.00026939906640532278, -0.00023779195323443418},
            {-0.0024466746444903481, 0.00022541051686156562, -0.00032192317838321434},
            {0.0019338240475487468, -8.8746368606900937e-05, 0.00049533492580329627},
            {-0.0018404229170903483, 7.5844093709869997e-05, -0.00049350901804910002},
            {0.0017156534860245892, 3.7345235972835756e-05, 0.0007411164716021292},
            {-0.0009460118308925611, 3.2811288545795183e-05, -0.00018460469050516324},
            {0.0014166799208792343, -0.00017080528579174628, -1.1234151255337116e-05}}};
    const auto geometry = fuelsim::elements::make_c3d8rt_geometry(coordinates);
    const bool native_increment_passed = test_reduced_finite_native_increment(geometry);
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = 300.0;
        for (std::size_t component = 0; component < 3; ++component)
            state[8 * (component + 1) + node] = displacement[node][component];
    }
    constexpr double young_modulus = 2.0e8, poisson_ratio = 0.25;
    const fuelsim::IsotropicThermoelasticMaterial material(
        fuelsim::test::thermoelastic(0.0, 1.0, young_modulus, poisson_ratio, 0.0, 300.0));
    fuelsim::CartesianMaterialHistory history(1);
    // A zero-increment hold with analytic isotropic compliance imposes the native
    // stress independently of this element's strain update. This checks force
    // assembly at the native state, not the preceding ten-step material history.
    const double trace = native_stress[0] + native_stress[1] + native_stress[2];
    for (std::size_t component = 0; component < 6; ++component)
        history[0].elastic_strain[component] =
            ((1.0 + poisson_ratio) * native_stress[component] - (component < 3 ? poisson_ratio * trace : 0.0))
            / young_modulus;
    const fuelsim::CartesianMaterialHistory zero_history(1);
    const fuelsim::elements::C3d8Input loaded_input{material,
        geometry,
        state,
        state,
        &history,
        1.0,
        1.0,
        0.0,
        fuelsim::StrainFormulation::finite,
        false,
        300.0};
    auto zero_input = loaded_input;
    zero_input.committed_history = &zero_history;
    const auto loaded = fuelsim::elements::evaluate_c3d8rt(loaded_input, {true, true, true, false});
    const auto passive = fuelsim::elements::evaluate_c3d8rt(loaded_input, {true, false, false, false});
    const auto baseline = fuelsim::elements::evaluate_c3d8rt(zero_input, {true, false, false, false});
    const auto& stress = loaded.history[0].stress;
    const std::array<double, 6> actual_stress = {stress.xx, stress.yy, stress.zz, stress.xy, stress.yz, stress.xz};
    bool passed = check(loaded.residual == passive.residual,
                      "native-state finite C3D8RT residual-only and Jacobian paths are identical")
                  && native_increment_passed;
    for (std::size_t component = 0; component < 6; ++component)
        passed = check(near(actual_stress[component], native_stress[component], 1.0e-13),
                     "analytic elastic compliance recovers every prescribed native stress component")
                 && passed;
    double difference_squared = 0.0, reference_squared = 0.0;
    double actual_peak = 0.0, reference_peak = 0.0, maximum_relative = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t row = 8 * (component + 1) + node;
            // The same-state zero-history subtraction cancels the production
            // hourglass force without calling its force or gradient helpers.
            const double actual = loaded.residual[row] - baseline.residual[row];
            const double reference = native_reaction[node][component] - weak_hourglass_force[node][component];
            const double difference = actual - reference;
            difference_squared += difference * difference;
            reference_squared += reference * reference;
            actual_peak = std::max(actual_peak, std::abs(actual));
            reference_peak = std::max(reference_peak, std::abs(reference));
            maximum_relative = std::max(maximum_relative, std::abs(difference / reference));
        }
    const double relative_l2 = std::sqrt(difference_squared / reference_squared);
    const double relative_peak = std::abs(actual_peak - reference_peak) / reference_peak;
    std::cout << "hex8_c3d8rt_native_bulk_force_relative_l2=" << relative_l2 << '\n'
              << "hex8_c3d8rt_native_bulk_force_relative_peak=" << relative_peak << '\n'
              << "hex8_c3d8rt_native_bulk_force_maximum_relative=" << maximum_relative << '\n';
    passed = check(std::max({relative_l2, relative_peak, maximum_relative}) < 1.0e-9,
                 "all 24 finite C3D8RT bulk-force components match the independent native-derived fixture")
             && passed;

    // Keep the same distorted reference element and nonaffine displacement mode,
    // now with a nonzero increment and nonuniform temperatures. E(T) and thermal
    // expansion activate the force-temperature chain; conductivity activates the
    // thermal-displacement chain together with the current-geometry capacity.
    auto previous = state;
    for (std::size_t node = 0; node < 8; ++node) {
        previous[node] = 305.0 + 0.4 * static_cast<double>(node);
        state[node] = 320.0 + 1.7 * static_cast<double>(node);
        for (std::size_t component = 0; component < 3; ++component)
            previous[8 * (component + 1) + node] *= 0.6;
    }
    const fuelsim::IsotropicThermoelasticMaterial varying_material(
        fuelsim::test::thermoelastic(500.0, 4.0, young_modulus, poisson_ratio, 1.0e-5, 300.0, -1.0e5));
    const fuelsim::elements::C3d8Input tangent_input{varying_material,
        geometry,
        state,
        previous,
        &history,
        2.0,
        7.0,
        3.0,
        fuelsim::StrainFormulation::finite,
        true,
        300.0};
    const auto tangent = fuelsim::elements::evaluate_c3d8rt(tangent_input, {true, true, false, false});
    const auto tangent_passive = fuelsim::elements::evaluate_c3d8rt(tangent_input, {true, false, false, false});
    passed = check(tangent.residual == tangent_passive.residual,
                 "distorted finite C3D8RT residual-only and Jacobian paths are identical")
             && passed;
    std::array<double, 4> block_errors{};
    std::array<double, 32> direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.31 * static_cast<double>(dof + 1));
    for (std::size_t block = 0; block < block_errors.size(); ++block) {
        const bool thermal_rows = block < 2, thermal_columns = block % 2 == 0;
        const double step = thermal_columns ? 1.0e-4 : 1.0e-7;
        auto plus = state, minus = state;
        for (std::size_t column = 0; column < direction.size(); ++column)
            if ((column < 8) == thermal_columns) {
                plus[column] += step * direction[column];
                minus[column] -= step * direction[column];
            }
        const auto evaluate_trial = [&](const fuelsim::Hex8LocalValues& trial) {
            const fuelsim::elements::C3d8Input input{tangent_input.material,
                tangent_input.geometry,
                trial,
                tangent_input.committed_state,
                tangent_input.committed_history,
                tangent_input.time_step,
                tangent_input.time,
                tangent_input.volumetric_heat_source,
                tangent_input.strain_formulation,
                tangent_input.include_thermal_time_term,
                tangent_input.initial_temperature};
            return fuelsim::elements::evaluate_c3d8rt(input, {true, false, false, false});
        };
        const auto plus_result = evaluate_trial(plus), minus_result = evaluate_trial(minus);
        double error_squared = 0.0, scale_squared = 0.0;
        for (std::size_t row = 0; row < 32; ++row) {
            if ((row < 8) != thermal_rows)
                continue;
            double analytic = 0.0;
            for (std::size_t column = 0; column < 32; ++column)
                if ((column < 8) == thermal_columns)
                    analytic += tangent.jacobian[row * 32 + column] * direction[column];
            const double numerical = (plus_result.residual[row] - minus_result.residual[row]) / (2.0 * step);
            error_squared += (analytic - numerical) * (analytic - numerical);
            scale_squared += numerical * numerical;
        }
        block_errors[block] = std::sqrt(error_squared / scale_squared);
        passed = check(scale_squared > 0.0 && block_errors[block] < 2.0e-6,
                     "each distorted nonaffine finite C3D8RT Jacobian block matches centered differences")
                 && passed;
    }
    std::cout << "hex8_c3d8rt_distorted_finite_ktt_relative_error=" << block_errors[0] << '\n'
              << "hex8_c3d8rt_distorted_finite_ktu_relative_error=" << block_errors[1] << '\n'
              << "hex8_c3d8rt_distorted_finite_kut_relative_error=" << block_errors[2] << '\n'
              << "hex8_c3d8rt_distorted_finite_kuu_relative_error=" << block_errors[3] << '\n';

    // Material particles retain their nonuniform nodal temperatures while the
    // nonaffine geometry changes. Every nodal temperature increment is zero, so
    // constant thermal expansion must add no strain despite changing volume weights.
    for (std::size_t node = 0; node < 8; ++node)
        previous[node] = state[node];
    std::array<fuelsim::elements::C3d8Result, 2> thermal_hold{};
    for (std::size_t variant = 0; variant < thermal_hold.size(); ++variant) {
        const double expansion = variant == 0 ? 0.0 : 1.0e-5;
        const fuelsim::IsotropicThermoelasticMaterial hold_material(
            fuelsim::test::thermoelastic(0.0, 1.0, young_modulus, poisson_ratio, expansion, 300.0));
        const fuelsim::elements::C3d8Input hold_input{hold_material,
            geometry,
            state,
            previous,
            &history,
            2.0,
            7.0,
            0.0,
            fuelsim::StrainFormulation::finite,
            false,
            300.0};
        thermal_hold[variant] = fuelsim::elements::evaluate_c3d8rt(hold_input, {true, false, true, false});
    }
    double hold_strain_error = 0.0, hold_stress_error = 0.0, hold_force_error = 0.0;
    const auto& unheated = thermal_hold[0].history[0];
    const auto& heated = thermal_hold[1].history[0];
    const std::array<double, 6> unheated_stress = {unheated.stress.xx,
        unheated.stress.yy,
        unheated.stress.zz,
        unheated.stress.xy,
        unheated.stress.yz,
        unheated.stress.xz};
    const std::array<double, 6> heated_stress =
        {heated.stress.xx, heated.stress.yy, heated.stress.zz, heated.stress.xy, heated.stress.yz, heated.stress.xz};
    for (std::size_t component = 0; component < 6; ++component) {
        hold_strain_error = std::max(hold_strain_error,
            std::abs(heated.elastic_strain[component] - unheated.elastic_strain[component]));
        hold_stress_error =
            std::max(hold_stress_error, std::abs(heated_stress[component] - unheated_stress[component]));
    }
    for (std::size_t row = 8; row < 32; ++row)
        hold_force_error =
            std::max(hold_force_error, std::abs(thermal_hold[1].residual[row] - thermal_hold[0].residual[row]));
    std::cout << "hex8_c3d8rt_unchanged_nodal_temperature_strain_error=" << hold_strain_error << '\n'
              << "hex8_c3d8rt_unchanged_nodal_temperature_stress_error=" << hold_stress_error << '\n'
              << "hex8_c3d8rt_unchanged_nodal_temperature_force_error=" << hold_force_error << '\n';
    passed = check(hold_strain_error < 1.0e-14 && hold_stress_error < 1.0e-7 && hold_force_error < 1.0e-7,
                 "nonaffine motion with unchanged nodal temperatures creates no thermal eigenstrain increment")
             && passed;
    return passed;
}

bool test_reduced_native_thermal_hourglass() {
    // Native fixtures: verification/abaqus/
    // b528_hex8_c3d8rt_warped_operator_nodal.csv, step BASE;
    // b532_hex8_c3d8rt_finite_warped_operator_nodal.csv, case BASE.
    // Both prescribe every temperature and displacement, with constant k = 4,
    // no thermal time term and no source. B532 takes one finite increment from
    // the reference state. These RFL values are native output, not a formula
    // evaluated with the production geometry or stabilization helpers.
    const fuelsim::Hex8Coordinates coordinates = {{{0.0, 0.0, 0.0},
        {1.2, 0.1, -0.05},
        {1.1, 1.0, 0.1},
        {-0.1, 0.9, 0.0},
        {0.05, -0.05, 1.0},
        {1.15, 0.0, 1.2},
        {1.0, 1.1, 1.1},
        {-0.05, 1.0, 0.9}}};
    const std::array<double, 8> temperature = {360.0, 410.0, 445.0, 385.0, 470.0, 430.0, 515.0, 455.0};
    const std::array<std::array<std::array<double, 3>, 8>, 2> displacement = {
        {{{{-1.209190383308867e-28, -1.915540125176229e-28, -1.658946874404279e-28},
             {0.00012, -2.0e-5, 3.0e-5},
             {0.00017, 9.000000000000001e-5, -4.0e-5},
             {-4.0e-5, 0.00011, 2.0e-5},
             {3.0e-5, -5.0e-5, 0.00014},
             {0.0001, 4.0e-5, 0.0001},
             {0.00023, 0.00016, 0.00018},
             {-8.000000000000001e-5, 6.999999999999999e-5, 0.00011}}},
            {{{4.986668050104441e-27, 6.754496808268046e-27, 5.673083331541935e-27},
                {0.048, -0.008, 0.012},
                {0.068, 0.036, -0.016},
                {-0.016, 0.04399999999999999, 0.008},
                {0.012, -0.02, 0.056},
                {0.04, 0.016, 0.04},
                {0.09200000000000001, 0.064, 0.07199999999999999},
                {-0.032, 0.028, 0.04399999999999999}}}}};
    const std::array<std::array<double, 8>, 2> native_reaction = {{{-141.0588035029169,
                                                                       -80.0705712640537,
                                                                       14.3189224961852,
                                                                       -84.27878878633828,
                                                                       51.75962860064374,
                                                                       9.096437254818547,
                                                                       147.6051743296219,
                                                                       82.62800087203955},
        {-148.5830822589203,
            -85.11862070779424,
            14.34897828494331,
            -87.78188141888214,
            54.8638643438567,
            8.635321596719107,
            152.5731453086042,
            91.06227485147329}}};
    const auto geometry = fuelsim::elements::make_c3d8rt_geometry(coordinates);
    const fuelsim::IsotropicThermoelasticMaterial material(
        fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 1.2e-5, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0));
    const fuelsim::CartesianMaterialHistory history(1);
    bool passed = true;
    for (std::size_t variant = 0; variant < displacement.size(); ++variant) {
        fuelsim::Hex8LocalValues state{}, committed{};
        for (std::size_t node = 0; node < 8; ++node) {
            state[node] = temperature[node];
            committed[node] = 300.0;
            for (std::size_t component = 0; component < 3; ++component)
                state[8 * (component + 1) + node] = displacement[variant][node][component];
        }
        const auto formulation = variant == 0 ? fuelsim::StrainFormulation::small : fuelsim::StrainFormulation::finite;
        const fuelsim::elements::C3d8Input
            input{material, geometry, state, committed, &history, 1.0, 1.0, 0.0, formulation, false, 300.0};
        const auto passive = fuelsim::elements::evaluate_c3d8rt(input, {true, false, false, false});
        const auto active = fuelsim::elements::evaluate_c3d8rt(input, {true, true, false, false});
        double maximum_error = 0.0;
        for (std::size_t node = 0; node < 8; ++node)
            maximum_error = std::max(maximum_error, std::abs(passive.residual[node] - native_reaction[variant][node]));
        std::cout << "hex8_c3d8rt_native_thermal_" << (variant == 0 ? "small" : "finite")
                  << "_maximum_absolute_error=" << maximum_error << '\n';
        passed = check(maximum_error < 1.0e-9,
                     "all eight distorted C3D8RT thermal reactions match the independent native fixture")
                 && check(passive.residual == active.residual,
                     "native thermal fixtures have identical residual-only and Jacobian-path residuals")
                 && passed;
    }
    return passed;
}

int run_c3d8rt_tests() {
    bool passed = test_history_geometry(true);
    passed = test_fixed_initial_mass_capacity(true) && passed;
    passed = test_reduced_integration_inelastic_jacobian() && passed;
    passed = test_reduced_integration_thermoelastic_capacity_gate() && passed;
    passed = test_reduced_integration_hourglass_energy() && passed;
    passed = test_reduced_finite_native_bulk_force() && passed;
    passed = test_reduced_native_thermal_hourglass() && passed;
    return passed ? 0 : 1;
}
} // namespace

int main() {
    return run_c3d8rt_tests();
}
