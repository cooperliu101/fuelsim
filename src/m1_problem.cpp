#include "fuelsim/m1_problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim {
namespace {

std::size_t checked_node_count(const StructuredRzMesh& fuel,
                               const StructuredRzMesh& cladding) {
    const std::size_t fuel_nodes = fuel.nodes().size();
    const std::size_t cladding_nodes = cladding.nodes().size();
    if (fuel_nodes > std::numeric_limits<std::size_t>::max() - cladding_nodes)
        throw std::length_error("M1Problem combined node count overflows");
    return fuel_nodes + cladding_nodes;
}

Quad4Coordinates element_coordinates(const StructuredRzMesh& mesh,
                                     const Quad4Element& element) {
    Quad4Coordinates coordinates{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node)
        coordinates[node] = mesh.nodes().at(element.nodes[node]);
    return coordinates;
}

std::array<std::size_t, 4> global_element_nodes(const Quad4Element& element,
                                                std::size_t node_offset) {
    std::array<std::size_t, 4> nodes{};
    for (std::size_t node = 0; node < nodes.size(); ++node)
        nodes[node] = node_offset + element.nodes[node];
    return nodes;
}

void append_boundary_conditions(const std::vector<std::size_t>& local_nodes,
                                std::size_t node_offset, Field field,
                                double value, const DofMap& dof_map,
                                std::vector<DirichletCondition>& output) {
    for (std::size_t local_node : local_nodes) {
        output.push_back({dof_map.dof(field, node_offset + local_node), value});
    }
}

void validate_dirichlet_conditions(
    std::vector<DirichletCondition>& conditions) {
    std::sort(conditions.begin(), conditions.end(),
              [](const DirichletCondition& lhs, const DirichletCondition& rhs) {
                  return lhs.dof < rhs.dof;
              });

    for (std::size_t index = 1; index < conditions.size(); ++index) {
        const DirichletCondition& previous = conditions[index - 1];
        const DirichletCondition& current = conditions[index];
        if (previous.dof != current.dof)
            continue;
        if (previous.value != current.value)
            throw std::invalid_argument(
                "M1Problem has conflicting Dirichlet conditions");
        throw std::invalid_argument(
            "M1Problem has duplicate Dirichlet conditions");
    }
}

} // namespace

M1Problem::M1Problem(M1Parameters parameters)
    : parameters_(parameters),
      fuel_mesh_(StructuredRzMesh::make_annulus(
          0.0, parameters.fuel_radius, parameters.length,
          parameters.fuel_radial_elements, parameters.axial_elements)),
      cladding_mesh_(StructuredRzMesh::make_annulus(
          parameters.cladding_inner_radius, parameters.cladding_outer_radius,
          parameters.length, parameters.cladding_radial_elements,
          parameters.axial_elements)),
      dof_map_(checked_node_count(fuel_mesh_, cladding_mesh_)),
      fuel_kernel_(IsotropicThermoelasticMaterial(parameters.fuel),
                   parameters.volumetric_heat_source),
      cladding_kernel_(IsotropicThermoelasticMaterial(parameters.cladding),
                       0.0),
      interface_kernel_({parameters.gap_conductivity, parameters.minimum_gap,
                         parameters.contact_penalty}) {
    if (!(parameters_.cladding_inner_radius > parameters_.fuel_radius))
        throw std::invalid_argument(
            "M1Problem requires a positive initial fuel-cladding gap");
    if (!std::isfinite(parameters_.volumetric_heat_source) ||
        !(parameters_.volumetric_heat_source >= 0.0))
        throw std::invalid_argument(
            "M1Problem volumetric_heat_source must be finite and "
            "nonnegative");
    if (!(parameters_.gap_conductivity > 0.0))
        throw std::invalid_argument(
            "M1Problem gap_conductivity must be positive");
    if (!(parameters_.contact_penalty > 0.0))
        throw std::invalid_argument(
            "M1Problem contact_penalty must be positive");
    if (!std::isfinite(parameters_.outer_temperature) ||
        !(parameters_.outer_temperature > 0.0) ||
        !std::isfinite(parameters_.initial_temperature) ||
        !(parameters_.initial_temperature > 0.0))
        throw std::invalid_argument(
            "M1Problem temperatures must be finite and positive");

    build_geometries();
    build_dirichlet_conditions();
}

const M1Parameters& M1Problem::parameters() const noexcept {
    return parameters_;
}

const StructuredRzMesh& M1Problem::fuel_mesh() const noexcept {
    return fuel_mesh_;
}

const StructuredRzMesh& M1Problem::cladding_mesh() const noexcept {
    return cladding_mesh_;
}

const DofMap& M1Problem::dof_map() const noexcept {
    return dof_map_;
}

const Quad4RzThermoelasticKernel& M1Problem::fuel_kernel() const noexcept {
    return fuel_kernel_;
}

const Quad4RzThermoelasticKernel& M1Problem::cladding_kernel() const noexcept {
    return cladding_kernel_;
}

const Line2RzGapContactKernel& M1Problem::interface_kernel() const noexcept {
    return interface_kernel_;
}

std::size_t M1Problem::fuel_node_count() const noexcept {
    return fuel_mesh_.nodes().size();
}

std::size_t M1Problem::cladding_node_offset() const noexcept {
    return fuel_node_count();
}

std::size_t M1Problem::fuel_global_node(std::size_t local_node) const {
    if (local_node >= fuel_node_count())
        throw std::out_of_range("M1Problem fuel node is out of range");
    return local_node;
}

std::size_t M1Problem::cladding_global_node(std::size_t local_node) const {
    if (local_node >= cladding_mesh_.nodes().size())
        throw std::out_of_range("M1Problem cladding node is out of range");
    return cladding_node_offset() + local_node;
}

std::size_t M1Problem::dof_count() const noexcept {
    return dof_map_.dof_count();
}

std::size_t M1Problem::contribution_count() const noexcept {
    return fuel_element_count() + cladding_element_count() + interface_count();
}

std::size_t M1Problem::fuel_element_count() const noexcept {
    return fuel_mesh_.elements().size();
}

std::size_t M1Problem::cladding_element_count() const noexcept {
    return cladding_mesh_.elements().size();
}

std::size_t M1Problem::interface_count() const noexcept {
    return interface_geometries_.size();
}

const std::vector<DirichletCondition>&
M1Problem::dirichlet_conditions() const noexcept {
    return dirichlet_conditions_;
}

std::vector<double> M1Problem::initial_state() const {
    std::vector<double> state(dof_count(), 0.0);
    for (std::size_t node = 0; node < dof_map_.node_count(); ++node)
        state[dof_map_.temperature(node)] = parameters_.initial_temperature;

    for (const DirichletCondition& condition : dirichlet_conditions_)
        state[condition.dof] = condition.value;
    return state;
}

LocalDofs M1Problem::contribution_dofs(std::size_t contribution_index) const {
    if (contribution_index < fuel_element_count())
        return fuel_element_dofs(contribution_index);

    contribution_index -= fuel_element_count();
    if (contribution_index < cladding_element_count())
        return cladding_element_dofs(contribution_index);

    contribution_index -= cladding_element_count();
    return interface_dofs(contribution_index);
}

LocalResidual M1Problem::contribution_residual(std::size_t contribution_index,
                                               const LocalValues& state) const {
    if (contribution_index < fuel_element_count()) {
        return fuel_kernel_.residual(fuel_geometries_.at(contribution_index),
                                     state);
    }

    contribution_index -= fuel_element_count();
    if (contribution_index < cladding_element_count()) {
        return cladding_kernel_.residual(
            cladding_geometries_.at(contribution_index), state);
    }

    contribution_index -= cladding_element_count();
    return interface_kernel_.residual(
        interface_geometries_.at(contribution_index), state);
}

LocalSystem M1Problem::linearize_contribution(std::size_t contribution_index,
                                              const LocalValues& state) const {
    if (contribution_index < fuel_element_count()) {
        return fuel_kernel_.linearize(fuel_geometries_.at(contribution_index),
                                      state);
    }

    contribution_index -= fuel_element_count();
    if (contribution_index < cladding_element_count()) {
        return cladding_kernel_.linearize(
            cladding_geometries_.at(contribution_index), state);
    }

    contribution_index -= cladding_element_count();
    return interface_kernel_.linearize(
        interface_geometries_.at(contribution_index), state);
}

LocalDofs M1Problem::fuel_element_dofs(std::size_t element_index) const {
    const Quad4Element& element = fuel_mesh_.elements().at(element_index);
    return dof_map_.local_dofs(global_element_nodes(element, 0));
}

LocalDofs M1Problem::cladding_element_dofs(std::size_t element_index) const {
    const Quad4Element& element = cladding_mesh_.elements().at(element_index);
    return dof_map_.local_dofs(
        global_element_nodes(element, cladding_node_offset()));
}

LocalDofs M1Problem::interface_dofs(std::size_t interface_index) const {
    const Line2BoundaryElement& fuel_edge =
        fuel_mesh_.boundary_elements(BoundaryId::radial_outer)
            .at(interface_index);
    const Line2BoundaryElement& cladding_edge =
        cladding_mesh_.boundary_elements(BoundaryId::radial_inner)
            .at(interface_index);
    const std::array<std::size_t, 4> nodes = {
        fuel_global_node(fuel_edge.nodes[0]),
        fuel_global_node(fuel_edge.nodes[1]),
        cladding_global_node(cladding_edge.nodes[0]),
        cladding_global_node(cladding_edge.nodes[1]),
    };
    return dof_map_.local_dofs(nodes);
}

const Line2RzInterfaceGeometry&
M1Problem::interface_geometry(std::size_t interface_index) const {
    return interface_geometries_.at(interface_index);
}

InterfaceSummary
M1Problem::summarize_interface(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "M1Problem interface summary state size mismatch");

    InterfaceSummary summary = {
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
    };
    const std::size_t first_interface =
        fuel_element_count() + cladding_element_count();

    for (std::size_t interface = 0; interface < interface_count();
         ++interface) {
        const LocalValues local_state =
            contribution_state(first_interface + interface, state);
        const InterfaceQuadratureValues values =
            interface_kernel_.quadrature_values(
                interface_geometries_[interface], local_state);
        const auto& points = interface_geometries_[interface].points;

        for (std::size_t q = 0; q < values.size(); ++q) {
            summary.minimum_gap = std::min(summary.minimum_gap, values[q].gap);
            summary.maximum_gap = std::max(summary.maximum_gap, values[q].gap);
            summary.maximum_contact_pressure =
                std::max(summary.maximum_contact_pressure, values[q].pressure);
            summary.total_heat_rate +=
                points[q].weighted_measure * values[q].heat_flux;
            summary.total_contact_force +=
                points[q].weighted_measure * values[q].pressure;
        }
    }

    if (interface_count() == 0)
        throw std::logic_error("M1Problem has no interface contributions");
    return summary;
}

void M1Problem::add_state_independent_residual(
    std::vector<double>& residual) const {
    (void)residual;
}

void M1Problem::build_geometries() {
    fuel_geometries_.reserve(fuel_mesh_.elements().size());
    for (const Quad4Element& element : fuel_mesh_.elements()) {
        fuel_geometries_.push_back(
            make_quad4_rz_geometry(element_coordinates(fuel_mesh_, element)));
    }

    cladding_geometries_.reserve(cladding_mesh_.elements().size());
    for (const Quad4Element& element : cladding_mesh_.elements()) {
        cladding_geometries_.push_back(make_quad4_rz_geometry(
            element_coordinates(cladding_mesh_, element)));
    }

    const auto& fuel_edges =
        fuel_mesh_.boundary_elements(BoundaryId::radial_outer);
    const auto& cladding_edges =
        cladding_mesh_.boundary_elements(BoundaryId::radial_inner);
    if (fuel_edges.size() != cladding_edges.size())
        throw std::invalid_argument(
            "M1Problem requires matching axial interface meshes");

    interface_geometries_.reserve(fuel_edges.size());
    for (std::size_t edge = 0; edge < fuel_edges.size(); ++edge) {
        Line2InterfaceSideCoordinates fuel_coordinates{};
        Line2InterfaceSideCoordinates cladding_coordinates{};
        for (std::size_t node = 0; node < 2; ++node) {
            fuel_coordinates[node] =
                fuel_mesh_.nodes().at(fuel_edges[edge].nodes[node]);
            cladding_coordinates[node] =
                cladding_mesh_.nodes().at(cladding_edges[edge].nodes[node]);
        }
        interface_geometries_.push_back(make_line2_rz_interface_geometry(
            fuel_coordinates, cladding_coordinates));
    }
}

void M1Problem::build_dirichlet_conditions() {
    append_boundary_conditions(
        cladding_mesh_.boundary_nodes(BoundaryId::radial_outer),
        cladding_node_offset(), Field::temperature,
        parameters_.outer_temperature, dof_map_, dirichlet_conditions_);
    append_boundary_conditions(
        fuel_mesh_.boundary_nodes(BoundaryId::radial_inner), 0,
        Field::radial_displacement, 0.0, dof_map_, dirichlet_conditions_);
    append_boundary_conditions(fuel_mesh_.boundary_nodes(BoundaryId::bottom), 0,
                               Field::axial_displacement, 0.0, dof_map_,
                               dirichlet_conditions_);
    append_boundary_conditions(
        cladding_mesh_.boundary_nodes(BoundaryId::bottom),
        cladding_node_offset(), Field::axial_displacement, 0.0, dof_map_,
        dirichlet_conditions_);

    validate_dirichlet_conditions(dirichlet_conditions_);
}

} // namespace fuelsim
