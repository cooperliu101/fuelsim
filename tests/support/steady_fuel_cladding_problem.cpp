#include "support/steady_fuel_cladding_problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

std::size_t checked_node_count(const StructuredRzMesh& fuel,
                               const StructuredRzMesh& cladding) {
    const std::size_t fuel_nodes = fuel.nodes().size();
    const std::size_t cladding_nodes = cladding.nodes().size();
    if (fuel_nodes > std::numeric_limits<std::size_t>::max() - cladding_nodes)
        throw std::length_error(
            "SteadyFuelCladdingProblem combined node count overflows");
    return fuel_nodes + cladding_nodes;
}

Quad4Coordinates element_coordinates(const StructuredRzMesh& mesh,
                                     const Quad4Element& element) {
    Quad4Coordinates coordinates{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node)
        coordinates[node] = mesh.nodes().at(element.nodes[node]);
    return coordinates;
}

Line2InterfaceSideCoordinates
edge_coordinates(const StructuredRzMesh& mesh,
                 const Line2BoundaryElement& edge) {
    return {{
        mesh.nodes().at(edge.nodes[0]),
        mesh.nodes().at(edge.nodes[1]),
    }};
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
    for (std::size_t local_node : local_nodes)
        output.push_back({dof_map.dof(field, node_offset + local_node), value});
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
            throw std::invalid_argument("SteadyFuelCladdingProblem has "
                                        "conflicting Dirichlet conditions");
        throw std::invalid_argument(
            "SteadyFuelCladdingProblem has duplicate Dirichlet conditions");
    }
}

std::size_t
find_containing_segment(double z, const StructuredRzMesh& mesh,
                        const std::vector<Line2BoundaryElement>& edges) {
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        const double lower = mesh.nodes().at(edges[edge].nodes[0]).z;
        const double upper = mesh.nodes().at(edges[edge].nodes[1]).z;
        const bool includes_upper = edge + 1 == edges.size();
        if (z >= lower && (z < upper || (includes_upper && z <= upper)))
            return edge;
    }
    throw std::invalid_argument("SteadyFuelCladdingProblem interface point is "
                                "outside the cladding surface");
}

} // namespace

SteadyFuelCladdingProblem::SteadyFuelCladdingProblem(
    SteadyFuelCladdingParameters parameters)
    : SteadyFuelCladdingProblem(
          parameters,
          StructuredRzMesh::make_annulus(
              0.0, parameters.fuel_radius, parameters.fuel_length,
              parameters.fuel_radial_elements, parameters.axial_elements),
          StructuredRzMesh::make_annulus(
              parameters.cladding_inner_radius,
              parameters.cladding_outer_radius, parameters.cladding_length,
              parameters.cladding_radial_elements, parameters.axial_elements)) {
}

SteadyFuelCladdingProblem::SteadyFuelCladdingProblem(
    SteadyFuelCladdingParameters parameters, StructuredRzMesh fuel_mesh,
    StructuredRzMesh cladding_mesh)
    : _parameters(parameters), _fuel_mesh(std::move(fuel_mesh)),
      _cladding_mesh(std::move(cladding_mesh)),
      _dof_map(checked_node_count(_fuel_mesh, _cladding_mesh)),
      _fuel_kernel(IsotropicThermoelasticMaterial(parameters.fuel),
                   parameters.volumetric_heat_source,
                   StrainFormulation::small),
      _cladding_kernel(IsotropicThermoelasticMaterial(parameters.cladding),
                       0.0, StrainFormulation::small),
      _gap_heat_kernel({parameters.gap_conductivity, parameters.minimum_gap}),
      _contact_kernel({parameters.contact_penalty}) {
    const auto same_geometry = [](double actual, double expected) {
        const double scale =
            std::max({1.0, std::abs(actual), std::abs(expected)});
        return std::abs(actual - expected) <= 1.0e-12 * scale;
    };
    if (_fuel_mesh.radial_elements() != _parameters.fuel_radial_elements ||
        _cladding_mesh.radial_elements() !=
            _parameters.cladding_radial_elements ||
        _fuel_mesh.axial_elements() != _parameters.axial_elements ||
        _cladding_mesh.axial_elements() != _parameters.axial_elements ||
        !same_geometry(_fuel_mesh.inner_radius(), 0.0) ||
        !same_geometry(_fuel_mesh.outer_radius(), _parameters.fuel_radius) ||
        !same_geometry(_fuel_mesh.length(), _parameters.fuel_length) ||
        !same_geometry(_cladding_mesh.inner_radius(),
                       _parameters.cladding_inner_radius) ||
        !same_geometry(_cladding_mesh.outer_radius(),
                       _parameters.cladding_outer_radius) ||
        !same_geometry(_cladding_mesh.length(), _parameters.cladding_length))
        throw std::invalid_argument(
            "SteadyFuelCladdingProblem imported meshes do not match "
            "SteadyFuelCladdingParameters");
    if (!(_parameters.cladding_inner_radius > _parameters.fuel_radius))
        throw std::invalid_argument("SteadyFuelCladdingProblem requires a "
                                    "positive initial fuel-cladding gap");
    if (!std::isfinite(_parameters.fuel_length) ||
        !(_parameters.fuel_length > 0.0) ||
        !std::isfinite(_parameters.cladding_length) ||
        !(_parameters.cladding_length >= _parameters.fuel_length))
        throw std::invalid_argument("SteadyFuelCladdingProblem cladding_length "
                                    "must be at least fuel_length");
    if (!std::isfinite(_parameters.volumetric_heat_source) ||
        !(_parameters.volumetric_heat_source >= 0.0))
        throw std::invalid_argument("SteadyFuelCladdingProblem "
                                    "volumetric_heat_source must be finite and "
                                    "nonnegative");
    if (!(_parameters.gap_conductivity > 0.0))
        throw std::invalid_argument(
            "SteadyFuelCladdingProblem gap_conductivity must be positive");
    if (!(_parameters.contact_penalty > 0.0))
        throw std::invalid_argument(
            "SteadyFuelCladdingProblem contact_penalty must be positive");
    if (!std::isfinite(_parameters.outer_temperature) ||
        !(_parameters.outer_temperature > 0.0) ||
        !std::isfinite(_parameters.initial_temperature) ||
        !(_parameters.initial_temperature > 0.0))
        throw std::invalid_argument("SteadyFuelCladdingProblem temperatures "
                                    "must be finite and positive");

    build_geometries();
    build_dirichlet_conditions();
}

const SteadyFuelCladdingParameters&
SteadyFuelCladdingProblem::parameters() const noexcept {
    return _parameters;
}

const StructuredRzMesh& SteadyFuelCladdingProblem::fuel_mesh() const noexcept {
    return _fuel_mesh;
}

const StructuredRzMesh&
SteadyFuelCladdingProblem::cladding_mesh() const noexcept {
    return _cladding_mesh;
}

const DofMap& SteadyFuelCladdingProblem::dof_map() const noexcept {
    return _dof_map;
}

const Quad4RzThermoelasticKernel&
SteadyFuelCladdingProblem::fuel_kernel() const noexcept {
    return _fuel_kernel;
}

const Quad4RzThermoelasticKernel&
SteadyFuelCladdingProblem::cladding_kernel() const noexcept {
    return _cladding_kernel;
}

const Line2RzGapHeatKernel&
SteadyFuelCladdingProblem::gap_heat_kernel() const noexcept {
    return _gap_heat_kernel;
}

const NodeToLineRzContactKernel&
SteadyFuelCladdingProblem::contact_kernel() const noexcept {
    return _contact_kernel;
}

void SteadyFuelCladdingProblem::set_volumetric_heat_source(
    double volumetric_heat_source) {
    _fuel_kernel.set_volumetric_heat_source(volumetric_heat_source);
    _parameters.volumetric_heat_source = volumetric_heat_source;
}

std::size_t SteadyFuelCladdingProblem::fuel_node_count() const noexcept {
    return _fuel_mesh.nodes().size();
}

std::size_t SteadyFuelCladdingProblem::cladding_node_offset() const noexcept {
    return fuel_node_count();
}

std::size_t
SteadyFuelCladdingProblem::fuel_global_node(std::size_t local_node) const {
    if (local_node >= fuel_node_count())
        throw std::out_of_range(
            "SteadyFuelCladdingProblem fuel node is out of range");
    return local_node;
}

std::size_t
SteadyFuelCladdingProblem::cladding_global_node(std::size_t local_node) const {
    if (local_node >= _cladding_mesh.nodes().size())
        throw std::out_of_range(
            "SteadyFuelCladdingProblem cladding node is out of range");
    return cladding_node_offset() + local_node;
}

std::size_t SteadyFuelCladdingProblem::dof_count() const noexcept {
    return _dof_map.dof_count();
}

std::size_t SteadyFuelCladdingProblem::contribution_count() const noexcept {
    return fuel_element_count() + cladding_element_count() +
           thermal_interface_count() + contact_contribution_count();
}

std::size_t SteadyFuelCladdingProblem::fuel_element_count() const noexcept {
    return _fuel_mesh.elements().size();
}

std::size_t SteadyFuelCladdingProblem::cladding_element_count() const noexcept {
    return _cladding_mesh.elements().size();
}

std::size_t
SteadyFuelCladdingProblem::thermal_interface_count() const noexcept {
    return _thermal_interface_geometries.size();
}

std::size_t
SteadyFuelCladdingProblem::contact_contribution_count() const noexcept {
    return _contact_geometries.size();
}

const std::vector<DirichletCondition>&
SteadyFuelCladdingProblem::dirichlet_conditions() const noexcept {
    return _dirichlet_conditions;
}

std::vector<double> SteadyFuelCladdingProblem::initial_state() const {
    std::vector<double> state(dof_count(), 0.0);
    for (std::size_t node = 0; node < _dof_map.node_count(); ++node)
        state[_dof_map.temperature(node)] = _parameters.initial_temperature;

    for (const DirichletCondition& condition : _dirichlet_conditions)
        state[condition.dof] = condition.value;
    return state;
}

LocalDofs SteadyFuelCladdingProblem::contribution_dofs(
    std::size_t contribution_index) const {
    if (contribution_index < fuel_element_count())
        return fuel_element_dofs(contribution_index);

    contribution_index -= fuel_element_count();
    if (contribution_index < cladding_element_count())
        return cladding_element_dofs(contribution_index);

    contribution_index -= cladding_element_count();
    if (contribution_index < thermal_interface_count())
        return thermal_interface_dofs(contribution_index);

    contribution_index -= thermal_interface_count();
    return contact_contribution_dofs(contribution_index);
}

LocalResidual SteadyFuelCladdingProblem::contribution_residual(
    std::size_t contribution_index, const LocalValues& state) const {
    if (contribution_index < fuel_element_count())
        return _fuel_kernel.residual(_fuel_geometries.at(contribution_index),
                                     state);

    contribution_index -= fuel_element_count();
    if (contribution_index < cladding_element_count())
        return _cladding_kernel.residual(
            _cladding_geometries.at(contribution_index), state);

    contribution_index -= cladding_element_count();
    if (contribution_index < thermal_interface_count())
        return _gap_heat_kernel.residual(
            _thermal_interface_geometries.at(contribution_index), state);

    contribution_index -= thermal_interface_count();
    return _contact_kernel.residual(_contact_geometries.at(contribution_index),
                                    state);
}

LocalSystem SteadyFuelCladdingProblem::linearize_contribution(
    std::size_t contribution_index, const LocalValues& state) const {
    if (contribution_index < fuel_element_count())
        return _fuel_kernel.linearize(_fuel_geometries.at(contribution_index),
                                      state);

    contribution_index -= fuel_element_count();
    if (contribution_index < cladding_element_count())
        return _cladding_kernel.linearize(
            _cladding_geometries.at(contribution_index), state);

    contribution_index -= cladding_element_count();
    if (contribution_index < thermal_interface_count())
        return _gap_heat_kernel.linearize(
            _thermal_interface_geometries.at(contribution_index), state);

    contribution_index -= thermal_interface_count();
    return _contact_kernel.linearize(_contact_geometries.at(contribution_index),
                                     state);
}

LocalDofs
SteadyFuelCladdingProblem::fuel_element_dofs(std::size_t element_index) const {
    const Quad4Element& element = _fuel_mesh.elements().at(element_index);
    return _dof_map.local_dofs(global_element_nodes(element, 0));
}

LocalDofs SteadyFuelCladdingProblem::cladding_element_dofs(
    std::size_t element_index) const {
    const Quad4Element& element = _cladding_mesh.elements().at(element_index);
    return _dof_map.local_dofs(
        global_element_nodes(element, cladding_node_offset()));
}

LocalDofs SteadyFuelCladdingProblem::thermal_interface_dofs(
    std::size_t interface_index) const {
    return _dof_map.local_dofs(_thermal_interface_nodes.at(interface_index));
}

LocalDofs SteadyFuelCladdingProblem::contact_contribution_dofs(
    std::size_t contact_index) const {
    return _dof_map.local_dofs(_contact_contribution_nodes.at(contact_index));
}

const Quad4RzGeometry& SteadyFuelCladdingProblem::fuel_element_geometry(
    std::size_t element_index) const {
    return _fuel_geometries.at(element_index);
}

const Quad4RzGeometry& SteadyFuelCladdingProblem::cladding_element_geometry(
    std::size_t element_index) const {
    return _cladding_geometries.at(element_index);
}

const Line2RzHeatGeometry&
SteadyFuelCladdingProblem::thermal_interface_geometry(
    std::size_t interface_index) const {
    return _thermal_interface_geometries.at(interface_index);
}

const NodeToLineRzContactGeometry&
SteadyFuelCladdingProblem::contact_contribution_geometry(
    std::size_t contact_index) const {
    return _contact_geometries.at(contact_index);
}

std::vector<ContactNodeSummary>
SteadyFuelCladdingProblem::summarize_contact_nodes(
    const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SteadyFuelCladdingProblem contact summary state size mismatch");

    const auto& fuel_nodes =
        _fuel_mesh.boundary_nodes(BoundaryId::radial_outer);
    std::vector<ContactNodeSummary> result;
    result.reserve(fuel_nodes.size());
    for (std::size_t node : fuel_nodes) {
        result.push_back({
            _fuel_mesh.nodes().at(node).z,
            false,
            std::numeric_limits<double>::infinity(),
            0.0,
            0.0,
            0.0,
            0.0,
        });
    }

    const std::size_t first_contact = fuel_element_count() +
                                      cladding_element_count() +
                                      thermal_interface_count();
    for (std::size_t contact = 0; contact < contact_contribution_count();
         ++contact) {
        const LocalValues local_state =
            contribution_state(first_contact + contact, state);
        const ContactPointValue value =
            _contact_kernel.value(_contact_geometries[contact], local_state);
        if (!value.projected)
            continue;

        ContactNodeSummary& node =
            result.at(_contact_secondary_axial_indices[contact]);
        node.projected = true;
        node.gap = std::min(node.gap, value.gap);
        node.tributary_area += value.tributary_area;
        node.tributary_length += value.tributary_length;
        node.contact_force += value.contact_force;
    }

    for (ContactNodeSummary& node : result) {
        if (node.tributary_area > 0.0)
            node.pressure = node.contact_force / node.tributary_area;
    }
    return result;
}

InterfaceSummary SteadyFuelCladdingProblem::summarize_interface(
    const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SteadyFuelCladdingProblem interface summary state size mismatch");

    InterfaceSummary summary = {
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
        0,
        0,
        0.0,
    };
    const std::size_t first_thermal =
        fuel_element_count() + cladding_element_count();
    for (std::size_t interface = 0; interface < thermal_interface_count();
         ++interface) {
        const LocalValues local_state =
            contribution_state(first_thermal + interface, state);
        const HeatQuadratureValues values = _gap_heat_kernel.quadrature_values(
            _thermal_interface_geometries[interface], local_state);
        for (const HeatQuadratureValue& value : values) {
            summary.minimum_gap = std::min(summary.minimum_gap, value.gap);
            summary.maximum_gap = std::max(summary.maximum_gap, value.gap);
            summary.total_heat_rate += value.weighted_measure * value.heat_flux;
        }
    }

    const std::vector<ContactNodeSummary> contact_nodes =
        summarize_contact_nodes(state);
    for (const ContactNodeSummary& node : contact_nodes) {
        if (!node.projected)
            continue;
        ++summary.projected_contact_nodes;
        summary.minimum_contact_gap =
            std::min(summary.minimum_contact_gap, node.gap);
        summary.maximum_contact_pressure =
            std::max(summary.maximum_contact_pressure, node.pressure);
        summary.total_contact_force += node.contact_force;
        if (node.pressure > 0.0) {
            ++summary.active_contact_nodes;
            summary.active_contact_length += node.tributary_length;
        }
    }

    if (thermal_interface_count() == 0 || contact_contribution_count() == 0)
        throw std::logic_error(
            "SteadyFuelCladdingProblem has no interface contributions");
    return summary;
}

void SteadyFuelCladdingProblem::add_state_independent_residual(
    std::vector<double>& residual) const {
    (void)residual;
}

void SteadyFuelCladdingProblem::build_geometries() {
    _fuel_geometries.reserve(_fuel_mesh.elements().size());
    for (const Quad4Element& element : _fuel_mesh.elements())
        _fuel_geometries.push_back(
            make_quad4_rz_geometry(element_coordinates(_fuel_mesh, element)));

    _cladding_geometries.reserve(_cladding_mesh.elements().size());
    for (const Quad4Element& element : _cladding_mesh.elements())
        _cladding_geometries.push_back(make_quad4_rz_geometry(
            element_coordinates(_cladding_mesh, element)));

    const auto& fuel_edges =
        _fuel_mesh.boundary_elements(BoundaryId::radial_outer);
    const auto& cladding_edges =
        _cladding_mesh.boundary_elements(BoundaryId::radial_inner);

    constexpr double gauss = 0.577350269189625764509148780501957456;
    _thermal_interface_nodes.reserve(fuel_edges.size());
    _thermal_interface_geometries.reserve(fuel_edges.size());
    for (const Line2BoundaryElement& fuel_edge : fuel_edges) {
        const Line2InterfaceSideCoordinates secondary_coordinates =
            edge_coordinates(_fuel_mesh, fuel_edge);
        const double lower_gauss_z =
            0.5 * (1.0 + gauss) * secondary_coordinates[0].z +
            0.5 * (1.0 - gauss) * secondary_coordinates[1].z;
        const double upper_gauss_z =
            0.5 * (1.0 - gauss) * secondary_coordinates[0].z +
            0.5 * (1.0 + gauss) * secondary_coordinates[1].z;
        const std::size_t lower_segment = find_containing_segment(
            lower_gauss_z, _cladding_mesh, cladding_edges);
        const std::size_t upper_segment = find_containing_segment(
            upper_gauss_z, _cladding_mesh, cladding_edges);
        if (lower_segment != upper_segment)
            throw std::invalid_argument("SteadyFuelCladdingProblem requires "
                                        "each fuel edge's thermal Gauss points "
                                        "to project to one cladding segment");

        const Line2BoundaryElement& cladding_edge =
            cladding_edges[lower_segment];
        const Line2InterfaceSideCoordinates primary_coordinates =
            edge_coordinates(_cladding_mesh, cladding_edge);
        _thermal_interface_nodes.push_back({
            fuel_global_node(fuel_edge.nodes[0]),
            fuel_global_node(fuel_edge.nodes[1]),
            cladding_global_node(cladding_edge.nodes[0]),
            cladding_global_node(cladding_edge.nodes[1]),
        });
        _thermal_interface_geometries.push_back(make_line2_rz_heat_geometry(
            secondary_coordinates, primary_coordinates));
    }

    for (std::size_t fuel_edge_index = 0; fuel_edge_index < fuel_edges.size();
         ++fuel_edge_index) {
        const Line2BoundaryElement& fuel_edge = fuel_edges[fuel_edge_index];
        const Line2InterfaceSideCoordinates secondary_coordinates =
            edge_coordinates(_fuel_mesh, fuel_edge);

        for (std::size_t secondary = 0;
             secondary < line2_interface_side_node_count; ++secondary) {
            const double secondary_z = secondary_coordinates[secondary].z;
            const std::size_t containing = find_containing_segment(
                secondary_z, _cladding_mesh, cladding_edges);
            const std::size_t first_candidate =
                containing == 0 ? 0 : containing - 1;
            const std::size_t last_candidate =
                std::min(containing + 1, cladding_edges.size() - 1);

            for (std::size_t candidate = first_candidate;
                 candidate <= last_candidate; ++candidate) {
                const Line2BoundaryElement& cladding_edge =
                    cladding_edges[candidate];
                const Line2InterfaceSideCoordinates primary_coordinates =
                    edge_coordinates(_cladding_mesh, cladding_edge);
                _contact_contribution_nodes.push_back({
                    fuel_global_node(fuel_edge.nodes[0]),
                    fuel_global_node(fuel_edge.nodes[1]),
                    cladding_global_node(cladding_edge.nodes[0]),
                    cladding_global_node(cladding_edge.nodes[1]),
                });
                _contact_geometries.push_back(
                    make_node_to_line_rz_contact_geometry(
                        secondary_coordinates, primary_coordinates, secondary,
                        candidate == 0,
                        candidate + 1 == cladding_edges.size()));
                _contact_secondary_axial_indices.push_back(fuel_edge_index +
                                                           secondary);
            }
        }
    }
}

void SteadyFuelCladdingProblem::build_dirichlet_conditions() {
    append_boundary_conditions(
        _cladding_mesh.boundary_nodes(BoundaryId::radial_outer),
        cladding_node_offset(), Field::temperature,
        _parameters.outer_temperature, _dof_map, _dirichlet_conditions);
    append_boundary_conditions(
        _fuel_mesh.boundary_nodes(BoundaryId::radial_inner), 0,
        Field::radial_displacement, 0.0, _dof_map, _dirichlet_conditions);
    append_boundary_conditions(_fuel_mesh.boundary_nodes(BoundaryId::bottom), 0,
                               Field::axial_displacement, 0.0, _dof_map,
                               _dirichlet_conditions);
    append_boundary_conditions(
        _cladding_mesh.boundary_nodes(BoundaryId::bottom),
        cladding_node_offset(), Field::axial_displacement, 0.0, _dof_map,
        _dirichlet_conditions);

    validate_dirichlet_conditions(_dirichlet_conditions);
}

} // namespace fuelsim
