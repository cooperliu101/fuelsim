#include "fuelsim/problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

Quad4Coordinates element_coordinates(const StructuredRzMesh& mesh,
                                     const Quad4Element& element) {
    Quad4Coordinates coordinates{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node)
        coordinates[node] = mesh.nodes().at(element.nodes[node]);
    return coordinates;
}

void append_boundary_conditions(const std::vector<std::size_t>& nodes,
                                Field field, double value,
                                const DofMap& dof_map,
                                std::vector<DirichletCondition>& output) {
    for (std::size_t node : nodes)
        output.push_back({dof_map.dof(field, node), value});
}

} // namespace

M0Problem::M0Problem(M0Parameters parameters)
    : parameters_(parameters),
      mesh_(StructuredRzMesh::make_annulus(
          parameters.inner_radius, parameters.outer_radius, parameters.length,
          parameters.radial_elements, parameters.axial_elements)),
      dof_map_(mesh_.nodes().size()),
      kernel_(IsotropicThermoelasticMaterial(parameters.fuel),
              parameters.volumetric_heat_source) {
    if (!(parameters_.inner_pressure >= 0.0) ||
        !(parameters_.outer_pressure >= 0.0))
        throw std::invalid_argument("M0Problem pressures must be nonnegative");

    geometries_.reserve(mesh_.elements().size());
    for (const Quad4Element& element : mesh_.elements())
        geometries_.push_back(
            make_quad4_rz_geometry(element_coordinates(mesh_, element)));

    build_dirichlet_conditions();
}

const M0Parameters& M0Problem::parameters() const noexcept {
    return parameters_;
}

const StructuredRzMesh& M0Problem::mesh() const noexcept {
    return mesh_;
}

const DofMap& M0Problem::dof_map() const noexcept {
    return dof_map_;
}

const Quad4RzThermoelasticKernel& M0Problem::kernel() const noexcept {
    return kernel_;
}

std::size_t M0Problem::dof_count() const noexcept {
    return dof_map_.dof_count();
}

std::size_t M0Problem::element_count() const noexcept {
    return mesh_.elements().size();
}

std::size_t M0Problem::contribution_count() const noexcept {
    return element_count();
}

const std::vector<DirichletCondition>&
M0Problem::dirichlet_conditions() const noexcept {
    return dirichlet_conditions_;
}

std::vector<double> M0Problem::initial_state() const {
    std::vector<double> state(dof_count(), 0.0);
    for (std::size_t node = 0; node < mesh_.nodes().size(); ++node)
        state[dof_map_.temperature(node)] = parameters_.initial_temperature;

    for (const DirichletCondition& condition : dirichlet_conditions_)
        state[condition.dof] = condition.value;
    return state;
}

LocalDofs M0Problem::element_dofs(std::size_t element_index) const {
    return dof_map_.element_dofs(mesh_.elements().at(element_index));
}

LocalValues
M0Problem::element_state(std::size_t element_index,
                         const std::vector<double>& global_state) const {
    return contribution_state(element_index, global_state);
}

LocalResidual M0Problem::element_residual(std::size_t element_index,
                                          const LocalValues& state) const {
    return kernel_.residual(geometries_.at(element_index), state);
}

LocalSystem M0Problem::linearize_element(std::size_t element_index,
                                         const LocalValues& state) const {
    return kernel_.linearize(geometries_.at(element_index), state);
}

LocalDofs M0Problem::contribution_dofs(std::size_t contribution_index) const {
    return element_dofs(contribution_index);
}

LocalResidual M0Problem::contribution_residual(std::size_t contribution_index,
                                               const LocalValues& state) const {
    return element_residual(contribution_index, state);
}

LocalSystem M0Problem::linearize_contribution(std::size_t contribution_index,
                                              const LocalValues& state) const {
    return linearize_element(contribution_index, state);
}

void M0Problem::add_state_independent_residual(
    std::vector<double>& residual) const {
    add_pressure_residual(residual);
}

const Quad4RzGeometry&
M0Problem::element_geometry(std::size_t element_index) const {
    return geometries_.at(element_index);
}

void M0Problem::build_dirichlet_conditions() {
    append_boundary_conditions(
        mesh_.boundary_nodes(BoundaryId::radial_outer), Field::temperature,
        parameters_.outer_temperature, dof_map_, dirichlet_conditions_);

    if (parameters_.inner_radius == 0.0) {
        append_boundary_conditions(
            mesh_.boundary_nodes(BoundaryId::radial_inner),
            Field::radial_displacement, 0.0, dof_map_, dirichlet_conditions_);
    }

    append_boundary_conditions(mesh_.boundary_nodes(BoundaryId::bottom),
                               Field::axial_displacement, 0.0, dof_map_,
                               dirichlet_conditions_);

    std::sort(dirichlet_conditions_.begin(), dirichlet_conditions_.end(),
              [](const DirichletCondition& lhs, const DirichletCondition& rhs) {
                  return lhs.dof < rhs.dof;
              });

    for (std::size_t index = 1; index < dirichlet_conditions_.size(); ++index) {
        const DirichletCondition& previous = dirichlet_conditions_[index - 1];
        const DirichletCondition& current = dirichlet_conditions_[index];
        if (previous.dof != current.dof)
            continue;
        if (previous.value != current.value)
            throw std::invalid_argument(
                "M0Problem has conflicting Dirichlet conditions");
        throw std::invalid_argument(
            "M0Problem has duplicate Dirichlet conditions");
    }
}

void M0Problem::add_pressure_residual(std::vector<double>& residual) const {
    const std::array<double, 2> gauss = {
        -0.577350269189625764509148780501957456,
        0.577350269189625764509148780501957456,
    };

    const auto add_boundary = [&](BoundaryId boundary, double pressure,
                                  double normal_r) {
        if (pressure == 0.0)
            return;

        for (const Line2BoundaryElement& edge :
             mesh_.boundary_elements(boundary)) {
            const RzPoint& first = mesh_.nodes().at(edge.nodes[0]);
            const RzPoint& second = mesh_.nodes().at(edge.nodes[1]);
            const double dr = second.r - first.r;
            const double dz = second.z - first.z;
            const double line_jacobian = 0.5 * std::sqrt(dr * dr + dz * dz);

            for (double xi : gauss) {
                const std::array<double, 2> shape = {
                    0.5 * (1.0 - xi),
                    0.5 * (1.0 + xi),
                };
                const double radius = shape[0] * first.r + shape[1] * second.r;
                const double measure = 2.0 * pi * radius * line_jacobian;

                for (std::size_t node = 0; node < 2; ++node) {
                    const std::size_t dof =
                        dof_map_.radial_displacement(edge.nodes[node]);
                    residual[dof] +=
                        measure * pressure * normal_r * shape[node];
                }
            }
        }
    };

    add_boundary(BoundaryId::radial_inner, parameters_.inner_pressure, -1.0);
    add_boundary(BoundaryId::radial_outer, parameters_.outer_pressure, 1.0);
}

} // namespace fuelsim
