#include "fuelsim/steady_problem.hpp"

#include "spatial_assembly.hpp"

#include <memory>
#include <utility>

namespace fuelsim {

SteadyProblem::SteadyProblem(SpatialDefinition definition,
                             const UnstructuredQuad4Mesh& source_mesh)
    : _spatial(
          std::make_unique<SpatialAssembly>(std::move(definition), source_mesh)) {
    _region_kernels.reserve(_spatial->region_count());
    for (std::size_t region_value = 0; region_value < _spatial->region_count();
         ++region_value) {
        const RegionDefinition& value = _spatial->region(region_value);
        _region_kernels.emplace_back(
            IsotropicThermoelasticMaterial(value.material),
            _spatial->region_heat_source(region_value),
            value.strain_formulation);
    }
}

SteadyProblem::~SteadyProblem() = default;

const SpatialDefinition& SteadyProblem::definition() const noexcept {
    return _spatial->definition();
}

const DofMap& SteadyProblem::dof_map() const noexcept {
    return _spatial->dof_map();
}

std::size_t SteadyProblem::region_count() const noexcept {
    return _spatial->region_count();
}

std::size_t SteadyProblem::region_index(const std::string& name) const {
    return _spatial->region_index(name);
}

const RegionDefinition& SteadyProblem::region(std::size_t index) const {
    return _spatial->region(index);
}

const RegionMesh& SteadyProblem::region_mesh(std::size_t index) const {
    return _spatial->region_mesh(index);
}

const Quad4RzThermoelasticKernel&
SteadyProblem::region_kernel(std::size_t index) const {
    return _region_kernels.at(index);
}

std::size_t SteadyProblem::region_node_offset(std::size_t index) const {
    return _spatial->region_node_offset(index);
}

std::size_t SteadyProblem::region_element_count(std::size_t index) const {
    return _spatial->region_element_count(index);
}

std::size_t SteadyProblem::region_element_offset(std::size_t index) const {
    return _spatial->region_element_offset(index);
}

std::size_t SteadyProblem::volume_contribution_count() const noexcept {
    return _spatial->volume_contribution_count();
}

SpatialContributionType
SteadyProblem::contribution_type(std::size_t contribution_index) const {
    return _spatial->contribution_type(contribution_index);
}

std::pair<std::size_t, std::size_t>
SteadyProblem::element_location(std::size_t contribution_index) const {
    return _spatial->element_location(contribution_index);
}

const Quad4RzGeometry& SteadyProblem::region_element_geometry(
    std::size_t region_value, std::size_t element_index) const {
    return _spatial->region_element_geometry(region_value, element_index);
}

std::size_t SteadyProblem::contact_count() const noexcept {
    return _spatial->contact_count();
}

const ContactDefinition& SteadyProblem::contact(std::size_t index) const {
    return _spatial->contact(index);
}

const std::vector<std::vector<ContactPointHistory>>&
SteadyProblem::committed_contact_histories() const noexcept {
    return _spatial->committed_contact_histories();
}

void SteadyProblem::commit_contact_state(const std::vector<double>& state) {
    _spatial->commit_contact_state(state);
}

bool SteadyProblem::uses_augmented_contact() const noexcept {
    return _spatial->uses_augmented_contact();
}

AugmentedContactUpdate SteadyProblem::update_augmented_contact_multipliers(
    const std::vector<double>& state, std::size_t completed_updates) {
    return _spatial->update_augmented_contact_multipliers(state,
                                                          completed_updates);
}

void SteadyProblem::restore_contact_state(
    const std::vector<double>& state,
    std::vector<std::vector<ContactPointHistory>> histories) {
    _spatial->restore_contact_state(state, std::move(histories));
}

void SteadyProblem::set_load_factor(double value) {
    _spatial->set_load_factor(value);
    refresh_region_heat_sources();
}

void SteadyProblem::refresh_region_heat_sources() {
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value)
        _region_kernels[region_value].set_volumetric_heat_source(
            _spatial->region_heat_source(region_value));
}

double SteadyProblem::load_factor() const noexcept {
    return _spatial->load_factor();
}

void SteadyProblem::set_time(double value) {
    _spatial->set_time(value);
    refresh_region_heat_sources();
}

std::vector<double> SteadyProblem::initial_state() const {
    return _spatial->initial_state();
}

std::vector<ContactNodeSummary> SteadyProblem::summarize_contact_nodes(
    std::size_t contact_index, const std::vector<double>& state) const {
    return _spatial->summarize_contact_nodes(contact_index, state);
}

std::vector<std::size_t> SteadyProblem::contact_secondary_source_nodes(
    std::size_t contact_index) const {
    return _spatial->contact_secondary_source_nodes(contact_index);
}

InterfaceSummary SteadyProblem::summarize_interface(
    std::size_t contact_index, const std::vector<double>& state) const {
    return _spatial->summarize_interface(contact_index, state);
}

std::size_t SteadyProblem::dof_count() const noexcept {
    return _spatial->dof_count();
}

std::size_t SteadyProblem::contribution_count() const noexcept {
    return _spatial->contribution_count();
}

const std::vector<DirichletCondition>&
SteadyProblem::dirichlet_conditions() const noexcept {
    return _spatial->dirichlet_conditions();
}

void SteadyProblem::validate_state(const std::vector<double>& state) const {
    _spatial->validate_state(state);
}

std::vector<std::size_t> SteadyProblem::required_state_dofs(
    std::size_t contribution_begin, std::size_t contribution_end) const {
    return _spatial->required_state_dofs(contribution_begin, contribution_end);
}

void SteadyProblem::validate_local_state(
    std::size_t contribution_begin, std::size_t contribution_end,
    const GlobalStateView& state) const {
    _spatial->validate_local_state(contribution_begin, contribution_end,
                                   state);
}

LocalDofs
SteadyProblem::contribution_dofs(std::size_t contribution_index) const {
    return _spatial->contribution_dofs(contribution_index);
}

LocalResidual SteadyProblem::contribution_residual(
    std::size_t contribution_index, const LocalValues& state) const {
    if (contribution_index < volume_contribution_count()) {
        const auto location = element_location(contribution_index);
        return _region_kernels[location.first].residual(
            region_element_geometry(location.first, location.second), state);
    }
    return _spatial->contribution_residual(contribution_index, state);
}

LocalSystem SteadyProblem::linearize_contribution(
    std::size_t contribution_index, const LocalValues& state) const {
    if (contribution_index < volume_contribution_count()) {
        const auto location = element_location(contribution_index);
        return _region_kernels[location.first].linearize(
            region_element_geometry(location.first, location.second), state);
    }
    return _spatial->linearize_contribution(contribution_index, state);
}

} // namespace fuelsim
