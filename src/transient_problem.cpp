#include "fuelsim/transient_problem.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

constexpr std::uint64_t fnv_offset = 14695981039346656037ULL;
constexpr std::uint64_t fnv_prime = 1099511628211ULL;

void hash_bytes(std::uint64_t& hash, const void* data, std::size_t size) {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= static_cast<std::uint64_t>(bytes[index]);
        hash *= fnv_prime;
    }
}

void hash_size(std::uint64_t& hash, std::size_t value) {
    const std::uint64_t encoded = static_cast<std::uint64_t>(value);
    hash_bytes(hash, &encoded, sizeof(encoded));
}

void hash_integer(std::uint64_t& hash, std::int64_t value) {
    hash_bytes(hash, &value, sizeof(value));
}

void hash_double(std::uint64_t& hash, double value) {
    std::uint64_t encoded = 0;
    static_assert(sizeof(encoded) == sizeof(value));
    std::memcpy(&encoded, &value, sizeof(value));
    hash_bytes(hash, &encoded, sizeof(encoded));
}

void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_size(hash, value.size());
    hash_bytes(hash, value.data(), value.size());
}

void hash_thermoelastic(std::uint64_t& hash,
                        const ThermoelasticProperties& material) {
    hash_double(hash, material.conductivity_inverse_temperature);
    hash_double(hash, material.conductivity_offset);
    hash_double(hash, material.young_modulus);
    hash_double(hash, material.poisson_ratio);
    hash_double(hash, material.thermal_expansion);
    hash_double(hash, material.reference_temperature);
    hash_double(hash, material.young_modulus_temperature_coefficient);
    hash_double(hash, material.poisson_ratio_temperature_coefficient);
    hash_double(hash, material.thermal_expansion_temperature_coefficient);
}

bool finite_stress(const AxisymmetricStressValues& stress) {
    return std::isfinite(stress.rr) && std::isfinite(stress.zz) &&
           std::isfinite(stress.hoop) && std::isfinite(stress.rz);
}

bool valid_material_state(const MaterialPointState& state) {
    for (std::size_t component = 0; component < 4; ++component) {
        if (!std::isfinite(state.elastic_strain[component]) ||
            !std::isfinite(state.plastic_strain[component]) ||
            !std::isfinite(state.creep_strain[component]))
            return false;
    }
    return std::isfinite(state.equivalent_plastic_strain) &&
           state.equivalent_plastic_strain >= 0.0 &&
           std::isfinite(state.equivalent_creep_strain) &&
           state.equivalent_creep_strain >= 0.0;
}

double stress_strain_inner_product(
    const AxisymmetricStressValues& stress,
    const std::array<double, 4>& strain) noexcept {
    return stress.rr * strain[0] + stress.zz * strain[1] +
           stress.hoop * strain[2] + 2.0 * stress.rz * strain[3];
}

std::array<double, 4> strain_difference(const std::array<double, 4>& current,
                                        const std::array<double, 4>& old) {
    std::array<double, 4> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        result[component] = current[component] - old[component];
    return result;
}

RegionInelasticSummary
summarize_history(const std::vector<Quad4MaterialHistory>& history) noexcept {
    RegionInelasticSummary summary{0.0, 0.0};
    for (const Quad4MaterialHistory& element : history) {
        for (const MaterialPointState& point : element) {
            summary.maximum_equivalent_plastic_strain =
                std::max(summary.maximum_equivalent_plastic_strain,
                         point.equivalent_plastic_strain);
            summary.maximum_equivalent_creep_strain =
                std::max(summary.maximum_equivalent_creep_strain,
                         point.equivalent_creep_strain);
        }
    }
    return summary;
}

void validate_definition(const TransientProblemDefinition& definition) {
    if (definition.regions.size() != definition.spatial.regions.size())
        throw std::invalid_argument(
            "TransientProblem requires one transient material per region");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        if (definition.regions[region].region !=
            definition.spatial.regions[region].name)
            throw std::invalid_argument(
                "TransientProblem region material order does not match the "
                "spatial regions");
    }
}

} // namespace

TransientProblem::TransientProblem(TransientProblemDefinition definition,
                                   const UnstructuredQuad4Mesh& source_mesh)
    : _definition(std::move(definition)),
      _spatial_model(_definition.spatial, source_mesh),
      _committed_solution(_spatial_model.initial_state()), _committed_time(0.0),
      _committed_load_factor(0.0), _active_time_step(0.0),
      _active_end_time(0.0), _active_load_factor(0.0),
      _time_step_active(false) {
    validate_definition(_definition);
    _region_kernels.reserve(region_count());
    _material_histories.resize(region_count());
    _material_stresses.resize(region_count());
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        _region_kernels.emplace_back(
            IsotropicInelasticMaterial(
                _definition.spatial.regions[region_value].material,
                _definition.regions[region_value].material),
            0.0,
            _definition.spatial.regions[region_value].strain_formulation);
        _material_histories[region_value].resize(
            _spatial_model.region_element_count(region_value));
        _material_stresses[region_value].resize(
            _spatial_model.region_element_count(region_value));
    }
    _spatial_model.set_load_factor(0.0);
    _committed_solution = _spatial_model.initial_state();
    _spatial_model.restore_contact_state(
        _committed_solution, _spatial_model.committed_contact_histories());
}

const TransientProblemDefinition&
TransientProblem::definition() const noexcept {
    return _definition;
}

const DofMap& TransientProblem::dof_map() const noexcept {
    return _spatial_model.dof_map();
}

std::size_t TransientProblem::region_count() const noexcept {
    return _definition.spatial.regions.size();
}

std::size_t TransientProblem::region_index(const std::string& name) const {
    return _spatial_model.region_index(name);
}

std::size_t
TransientProblem::region_node_offset(std::size_t region_value) const {
    return _spatial_model.region_node_offset(region_value);
}

const RegionDefinition& TransientProblem::region(std::size_t index) const {
    return _spatial_model.region(index);
}

const RegionMesh& TransientProblem::region_mesh(std::size_t index) const {
    return _spatial_model.region_mesh(index);
}

const Quad4RzTransientKernel&
TransientProblem::region_kernel(std::size_t index) const {
    return _region_kernels.at(index);
}

const Quad4RzGeometry&
TransientProblem::region_element_geometry(std::size_t region_value,
                                          std::size_t element_index) const {
    return _spatial_model.region_element_geometry(region_value, element_index);
}

const std::vector<double>&
TransientProblem::committed_solution() const noexcept {
    return _committed_solution;
}

double TransientProblem::committed_time() const noexcept {
    return _committed_time;
}

double TransientProblem::committed_load_factor() const noexcept {
    return _committed_load_factor;
}

bool TransientProblem::time_step_active() const noexcept {
    return _time_step_active;
}

double TransientProblem::active_time_step() const {
    require_active_time_step();
    return _active_time_step;
}

double TransientProblem::active_end_time() const {
    require_active_time_step();
    return _active_end_time;
}

std::uint64_t TransientProblem::committed_state_signature() const {
    std::uint64_t hash = fnv_offset;
    hash_size(hash, dof_count());
    hash_size(hash, region_count());
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const RegionDefinition& spatial =
            _definition.spatial.regions[region_value];
        const TransientInelasticProperties& transient =
            _definition.regions[region_value].material;
        hash_string(hash, spatial.name);
        hash_string(hash, spatial.block);
        hash_integer(hash, spatial.block_id);
        hash_thermoelastic(hash, spatial.material);
        hash_double(hash, spatial.volumetric_heat_source);
        hash_double(hash, spatial.initial_temperature);
        hash_string(hash, spatial.heat_source_function);
        hash_integer(hash,
                     static_cast<std::int64_t>(spatial.strain_formulation));
        hash_double(hash, transient.density);
        hash_double(hash, transient.specific_heat);
        hash_integer(hash, static_cast<std::int64_t>(transient.behavior));
        hash_double(hash, transient.creep.coefficient);
        hash_double(hash, transient.creep.reference_stress);
        hash_double(hash, transient.creep.stress_exponent);
        hash_double(hash,
                    transient.creep.coefficient_temperature_coefficient);
        hash_double(
            hash,
            transient.creep.reference_stress_temperature_coefficient);
        hash_double(hash,
                    transient.creep.stress_exponent_temperature_coefficient);
        hash_double(hash, transient.plasticity.yield_stress);
        hash_double(hash, transient.plasticity.isotropic_hardening_modulus);
        hash_double(
            hash,
            transient.plasticity.yield_stress_temperature_coefficient);
        hash_double(hash,
                    transient.plasticity.hardening_temperature_coefficient);

        const RegionMesh& mesh = region_mesh(region_value);
        hash_size(hash, mesh.nodes().size());
        for (const RzPoint& point : mesh.nodes()) {
            hash_double(hash, point.r);
            hash_double(hash, point.z);
        }
        hash_size(hash, mesh.elements().size());
        for (const Quad4Element& element : mesh.elements()) {
            for (const std::size_t node : element.nodes)
                hash_size(hash, node);
        }
        for (const std::size_t source : mesh.source_node_ids())
            hash_size(hash, source);
        for (const std::size_t source : mesh.source_element_ids())
            hash_size(hash, source);
    }
    for (const ContactDefinition& contact : _definition.spatial.contacts) {
        hash_string(hash, contact.name);
        hash_string(hash, contact.primary);
        hash_string(hash, contact.secondary);
        hash_integer(hash, contact.thermal ? 1 : 0);
        hash_integer(hash, contact.mechanical ? 1 : 0);
        hash_double(hash, contact.gap_conductivity);
        hash_double(hash, contact.minimum_gap);
        hash_double(hash, contact.penalty);
        hash_double(hash, contact.friction_coefficient);
    }
    for (const BoundaryConditionDefinition& boundary :
         _definition.spatial.boundary_conditions) {
        hash_string(hash, boundary.name);
        hash_integer(hash, static_cast<std::int64_t>(boundary.type));
        hash_string(hash, boundary.boundary);
        hash_integer(hash, static_cast<std::int64_t>(boundary.field));
        hash_double(hash, boundary.value);
        hash_integer(hash, boundary.scale_with_load ? 1 : 0);
        hash_string(hash, boundary.function);
        hash_double(hash, boundary.heat_transfer_coefficient);
        hash_double(hash, boundary.ambient_temperature);
        hash_string(hash, boundary.coefficient_function);
        hash_string(hash, boundary.ambient_temperature_function);
        hash_integer(hash, boundary.use_displaced_geometry ? 1 : 0);
    }
    for (const PiecewiseLinearTimeTable& table :
         _definition.spatial.time_tables) {
        hash_string(hash, table.name());
        hash_size(hash, table.times().size());
        for (std::size_t entry = 0; entry < table.times().size(); ++entry) {
            hash_double(hash, table.times()[entry]);
            hash_double(hash, table.values()[entry]);
        }
    }
    hash_size(hash, contribution_count());
    for (std::size_t contribution = 0; contribution < contribution_count();
         ++contribution) {
        const LocalDofs dofs = contribution_dofs(contribution);
        for (const std::size_t dof : dofs)
            hash_size(hash, dof);
    }
    return hash;
}

std::vector<double> TransientProblem::time_events() const {
    std::vector<double> result;
    for (const PiecewiseLinearTimeTable& table :
         _definition.spatial.time_tables)
        result.insert(result.end(), table.times().begin(), table.times().end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

TransientCommittedState TransientProblem::committed_state() const {
    return {_committed_solution, _material_histories, _material_stresses,
            _spatial_model.committed_contact_histories(),
            _last_conservation_summary, _committed_time,
            _committed_load_factor};
}

void TransientProblem::restore_committed_state(TransientCommittedState state) {
    if (_time_step_active)
        throw std::logic_error(
            "TransientProblem cannot restore during an active time step");
    if (state.solution.size() != dof_count() ||
        state.material_histories.size() != region_count() ||
        state.material_stresses.size() != region_count() ||
        state.contact_histories.size() != _definition.spatial.contacts.size())
        throw std::invalid_argument(
            "Transient committed state layout does not match the problem");
    if (!std::isfinite(state.time) || state.time < 0.0 ||
        !std::isfinite(state.load_factor) || state.load_factor < 0.0)
        throw std::invalid_argument(
            "Transient committed time and load factor must be valid");
    for (std::size_t node = 0; node < dof_map().node_count(); ++node) {
        const double temperature =
            state.solution.at(dof_map().temperature(node));
        if (!std::isfinite(temperature) || !(temperature > 0.0) ||
            !std::isfinite(
                state.solution.at(dof_map().radial_displacement(node))) ||
            !std::isfinite(
                state.solution.at(dof_map().axial_displacement(node))))
            throw std::invalid_argument(
                "Transient committed nodal state must be finite with "
                "positive temperatures");
    }
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const std::size_t elements =
            region_mesh(region_value).elements().size();
        if (state.material_histories[region_value].size() != elements ||
            state.material_stresses[region_value].size() != elements)
            throw std::invalid_argument(
                "Transient committed element state layout does not match");
        for (std::size_t element = 0; element < elements; ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                if (!valid_material_state(
                        state.material_histories[region_value][element][q]) ||
                    !finite_stress(
                        state.material_stresses[region_value][element][q]))
                    throw std::invalid_argument(
                        "Transient committed integration-point state is "
                        "invalid");
            }
        }
    }

    _spatial_model.restore_contact_state(state.solution,
                                         state.contact_histories);
    _committed_solution = std::move(state.solution);
    _material_histories = std::move(state.material_histories);
    _material_stresses = std::move(state.material_stresses);
    _last_conservation_summary = state.conservation;
    _committed_time = state.time;
    _committed_load_factor = state.load_factor;
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _active_load_factor = _committed_load_factor;
    _spatial_model.set_time(_committed_time);
    _spatial_model.set_load_factor(_committed_load_factor);
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value)
        _region_kernels[region_value].set_volumetric_heat_source(
            _spatial_model.region_kernel(region_value)
                .volumetric_heat_source());
}

void TransientProblem::begin_time_step(const TransientStepInput& input) {
    if (_time_step_active)
        throw std::logic_error(
            "TransientProblem already has an active time step");
    if (!std::isfinite(input.end_time) || input.end_time <= _committed_time)
        throw std::invalid_argument(
            "TransientProblem end time must exceed committed time");
    if (!std::isfinite(input.load_factor) || input.load_factor < 0.0)
        throw std::invalid_argument(
            "TransientProblem load factor must be finite and nonnegative");
    _active_time_step = input.end_time - _committed_time;
    _active_end_time = input.end_time;
    _active_load_factor = input.load_factor;
    try {
        _spatial_model.set_time(input.end_time);
        _spatial_model.set_load_factor(input.load_factor);
    } catch (...) {
        _spatial_model.set_time(_committed_time);
        _spatial_model.set_load_factor(_committed_load_factor);
        _active_time_step = 0.0;
        _active_end_time = _committed_time;
        _active_load_factor = _committed_load_factor;
        throw;
    }
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value)
        _region_kernels[region_value].set_volumetric_heat_source(
            _spatial_model.region_kernel(region_value)
                .volumetric_heat_source());
    _time_step_active = true;
}

void TransientProblem::commit_time_step(
    const std::vector<double>& converged_solution) {
    require_active_time_step();
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument(
            "TransientProblem committed solution size mismatch");
    if (!std::all_of(converged_solution.begin(), converged_solution.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::domain_error(
            "TransientProblem committed solution must be finite");
    for (std::size_t node = 0; node < dof_map().node_count(); ++node) {
        if (!(converged_solution[dof_map().temperature(node)] > 0.0))
            throw std::domain_error(
                "TransientProblem committed temperatures must be positive");
    }

    std::vector<std::vector<Quad4MaterialHistory>> staged(region_count());
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>
        staged_stresses(region_count());
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        staged[region_value].resize(
            _spatial_model.region_element_count(region_value));
        staged_stresses[region_value].resize(
            _spatial_model.region_element_count(region_value));
        const std::size_t offset =
            _spatial_model.region_element_offset(region_value);
        for (std::size_t element = 0; element < staged[region_value].size();
             ++element) {
            const LocalValues state =
                contribution_state(offset + element, converged_solution);
            const LocalValues committed_state =
                contribution_state(offset + element, _committed_solution);
            staged[region_value][element] =
                _region_kernels[region_value].trial_state_values(
                    region_element_geometry(region_value, element), state,
                    committed_state,
                    _material_histories[region_value][element],
                    _active_time_step);
            staged_stresses[region_value][element] =
                _region_kernels[region_value].stress_values(
                    region_element_geometry(region_value, element), state,
                    committed_state,
                    _material_histories[region_value][element],
                    _active_time_step);
        }
    }

    const TransientConservationSummary conservation = summarize_active_step(
        converged_solution, staged, staged_stresses);
    _spatial_model.commit_contact_state(converged_solution);
    _material_histories.swap(staged);
    _material_stresses.swap(staged_stresses);
    _last_conservation_summary = conservation;
    _committed_solution = converged_solution;
    _committed_time = _active_end_time;
    _committed_load_factor = _active_load_factor;
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _time_step_active = false;
}

void TransientProblem::rollback_time_step() noexcept {
    if (!_time_step_active)
        return;
    _spatial_model.set_time(_committed_time);
    _spatial_model.set_load_factor(_committed_load_factor);
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value)
        _region_kernels[region_value].set_volumetric_heat_source(
            _spatial_model.region_kernel(region_value)
                .volumetric_heat_source());
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _active_load_factor = _committed_load_factor;
    _time_step_active = false;
}

const Quad4MaterialHistory&
TransientProblem::material_history(std::size_t region_value,
                                   std::size_t element_index) const {
    return _material_histories.at(region_value).at(element_index);
}

const std::array<AxisymmetricStressValues, 4>&
TransientProblem::material_stress(std::size_t region_value,
                                  std::size_t element_index) const {
    return _material_stresses.at(region_value).at(element_index);
}

RegionInelasticSummary
TransientProblem::summarize_region_history(std::size_t region_value) const {
    return summarize_history(_material_histories.at(region_value));
}

const TransientConservationSummary&
TransientProblem::last_conservation_summary() const noexcept {
    return _last_conservation_summary;
}

TransientConservationSummary TransientProblem::summarize_active_step(
    const std::vector<double>& converged_solution,
    const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
    const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>&
        staged_stresses) const {
    require_active_time_step();
    TransientConservationSummary result;
    std::vector<double> raw_residual(dof_count(), 0.0);

    for (std::size_t contribution = 0; contribution < contribution_count();
         ++contribution) {
        const LocalDofs dofs = contribution_dofs(contribution);
        const LocalValues state =
            contribution_state(contribution, converged_solution);
        const LocalResidual residual =
            contribution_residual(contribution, state);
        const SpatialContributionType type =
            _spatial_model.contribution_type(contribution);
        for (std::size_t local = 0; local < local_dof_count; ++local) {
            raw_residual[dofs[local]] += residual[local];
            const double increment =
                converged_solution[dofs[local]] -
                _committed_solution[dofs[local]];
            if (local < 4) {
                if (type == SpatialContributionType::thermal_contact)
                    result.interface_heat_imbalance += residual[local];
                else if (type == SpatialContributionType::convection)
                    result.convection_heat_rate += residual[local];
                continue;
            }
            const double work = residual[local] * increment;
            if (type == SpatialContributionType::volume)
                result.internal_mechanical_work_increment += work;
            else if (type == SpatialContributionType::mechanical_contact)
                result.contact_work_increment += work;
            else if (type == SpatialContributionType::pressure ||
                     type == SpatialContributionType::traction)
                result.pressure_traction_work_increment -= work;
        }
    }
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const Quad4RzTransientKernel& kernel =
            _region_kernels[region_value];
        const TransientInelasticProperties& properties = kernel.properties();
        const double heat_capacity = properties.density * properties.specific_heat;
        const std::size_t offset =
            _spatial_model.region_element_offset(region_value);
        for (std::size_t element = 0;
             element < staged_histories[region_value].size(); ++element) {
            const LocalValues current =
                contribution_state(offset + element, converged_solution);
            const LocalValues old =
                contribution_state(offset + element, _committed_solution);
            const Quad4RzGeometry& geometry =
                region_element_geometry(region_value, element);
            for (std::size_t q = 0; q < geometry.points.size(); ++q) {
                const RzQuadraturePoint& point = geometry.points[q];
                double current_temperature = 0.0;
                double old_temperature = 0.0;
                for (std::size_t node = 0; node < quad4_node_count; ++node) {
                    current_temperature += point.shape[node] * current[node];
                    old_temperature += point.shape[node] * old[node];
                }
                result.stored_heat_rate +=
                    point.weighted_measure * heat_capacity *
                    (current_temperature - old_temperature) /
                    _active_time_step;
                result.generated_heat_rate +=
                    point.weighted_measure * kernel.volumetric_heat_source();

                const MaterialPointState& old_history =
                    _material_histories[region_value][element][q];
                const MaterialPointState& new_history =
                    staged_histories[region_value][element][q];
                const AxisymmetricStressValues& old_stress =
                    _material_stresses[region_value][element][q];
                const AxisymmetricStressValues& new_stress =
                    staged_stresses[region_value][element][q];
                result.elastic_energy_change +=
                    0.5 * point.weighted_measure *
                    (stress_strain_inner_product(
                         new_stress, new_history.elastic_strain) -
                     stress_strain_inner_product(
                         old_stress, old_history.elastic_strain));
                result.plastic_dissipation_increment +=
                    point.weighted_measure * stress_strain_inner_product(
                        new_stress,
                        strain_difference(new_history.plastic_strain,
                                          old_history.plastic_strain));
                result.creep_dissipation_increment +=
                    point.weighted_measure * stress_strain_inner_product(
                        new_stress,
                        strain_difference(new_history.creep_strain,
                                          old_history.creep_strain));
            }
        }
    }

    std::vector<bool> constrained(dof_count(), false);
    for (const DirichletCondition& condition : dirichlet_conditions()) {
        constrained[condition.dof] = true;
        const double increment = converged_solution[condition.dof] -
                                 _committed_solution[condition.dof];
        if (condition.dof < dof_map().node_count())
            result.dirichlet_heat_input_rate += raw_residual[condition.dof];
        else
            result.dirichlet_reaction_work_increment +=
                raw_residual[condition.dof] * increment;
    }
    double thermal_residual_squared = 0.0;
    double mechanical_residual_squared = 0.0;
    for (std::size_t dof = 0; dof < dof_count(); ++dof) {
        if (constrained[dof])
            continue;
        const double value = raw_residual[dof];
        if (dof < dof_map().node_count())
            thermal_residual_squared += value * value;
        else
            mechanical_residual_squared += value * value;
    }
    result.unconstrained_thermal_residual_l2 =
        std::sqrt(thermal_residual_squared);
    result.unconstrained_mechanical_residual_l2 =
        std::sqrt(mechanical_residual_squared);
    result.global_thermal_balance =
        result.stored_heat_rate + result.convection_heat_rate +
        result.interface_heat_imbalance - result.generated_heat_rate -
        result.dirichlet_heat_input_rate;
    const double thermal_scale =
        std::abs(result.stored_heat_rate) +
        std::abs(result.convection_heat_rate) +
        std::abs(result.interface_heat_imbalance) +
        std::abs(result.generated_heat_rate) +
        std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0
            ? std::abs(result.global_thermal_balance) / thermal_scale
            : 0.0;
    result.mechanical_work_balance =
        result.internal_mechanical_work_increment +
        result.contact_work_increment -
        result.pressure_traction_work_increment -
        result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) +
        std::abs(result.contact_work_increment) +
        std::abs(result.pressure_traction_work_increment) +
        std::abs(result.dirichlet_reaction_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0
            ? std::abs(result.mechanical_work_balance) / mechanical_scale
            : 0.0;
    return result;
}

InterfaceSummary
TransientProblem::summarize_interface(std::size_t contact_index,
                                      const std::vector<double>& state) const {
    return _spatial_model.summarize_interface(contact_index, state);
}

std::vector<ContactNodeSummary> TransientProblem::summarize_contact_nodes(
    std::size_t contact_index, const std::vector<double>& state) const {
    return _spatial_model.summarize_contact_nodes(contact_index, state);
}

std::vector<std::size_t> TransientProblem::contact_secondary_source_nodes(
    std::size_t contact_index) const {
    return _spatial_model.contact_secondary_source_nodes(contact_index);
}

std::size_t TransientProblem::dof_count() const noexcept {
    return _spatial_model.dof_count();
}

std::size_t TransientProblem::contribution_count() const noexcept {
    return _spatial_model.contribution_count();
}

const std::vector<DirichletCondition>&
TransientProblem::dirichlet_conditions() const noexcept {
    return _spatial_model.dirichlet_conditions();
}

void TransientProblem::validate_state(const std::vector<double>& state) const {
    require_active_time_step();
    _spatial_model.validate_state(state);
}

LocalDofs
TransientProblem::contribution_dofs(std::size_t contribution_index) const {
    return _spatial_model.contribution_dofs(contribution_index);
}

LocalResidual
TransientProblem::contribution_residual(std::size_t contribution_index,
                                        const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _spatial_model.volume_contribution_count()) {
        const auto location =
            _spatial_model.element_location(contribution_index);
        return _region_kernels[location.first].residual(
            region_element_geometry(location.first, location.second), state,
            committed_element_state(contribution_index),
            _material_histories[location.first][location.second],
            _active_time_step);
    }
    return _spatial_model.contribution_residual(contribution_index, state);
}

LocalSystem
TransientProblem::linearize_contribution(std::size_t contribution_index,
                                         const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _spatial_model.volume_contribution_count()) {
        const auto location =
            _spatial_model.element_location(contribution_index);
        return _region_kernels[location.first].linearize(
            region_element_geometry(location.first, location.second), state,
            committed_element_state(contribution_index),
            _material_histories[location.first][location.second],
            _active_time_step);
    }
    return _spatial_model.linearize_contribution(contribution_index, state);
}

void TransientProblem::add_state_independent_residual(
    std::vector<double>& residual) const {
    (void)residual;
}

LocalValues TransientProblem::committed_element_state(
    std::size_t contribution_index) const {
    return contribution_state(contribution_index, _committed_solution);
}

void TransientProblem::require_active_time_step() const {
    if (!_time_step_active)
        throw std::logic_error("TransientProblem residual evaluation requires "
                               "an active time step");
}

} // namespace fuelsim
