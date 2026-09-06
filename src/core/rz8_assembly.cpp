#include "rz8_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim::rz8 {
SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad8Mesh& source)
    : SpatialLayout(definition, spatial_detail::resolve_block_ids(definition, source, true, true),
          spatial_detail::DofLayout::axisymmetric_rz) {
    std::vector<std::size_t> nodes, elements;
    std::vector<std::vector<std::size_t>> source_nodes;
    std::vector<std::vector<bool>> temperature_nodes;
    for (std::size_t r = 0; r < region_count(); ++r) {
        if (region(r).rz_element_formulation != RzElementFormulation::cax8t)
            throw std::invalid_argument("QUAD8 axisymmetric regions require element = cax8t");
        _meshes.push_back(Quad8RegionMesh::from_unstructured_block(source, _block_ids[r]));
        const auto& mesh = _meshes.back();
        nodes.push_back(mesh.nodes().size());
        elements.push_back(mesh.elements().size());
        source_nodes.push_back(mesh.source_node_ids());
        temperature_nodes.push_back(mesh.temperature_nodes());
        _geometries.emplace_back();
        for (const auto& e : mesh.elements()) {
            Quad8RzCoordinates coordinates;
            for (std::size_t n = 0; n < 8; ++n) coordinates[n] = mesh.nodes()[e.nodes[n]];
            _geometries.back().push_back(make_quad8_rz_geometry(coordinates));
        }
    }
    initialize_counts(nodes, elements);
    initialize_mixed_shared_nodes(source_nodes, temperature_nodes);
    for (std::size_t b = 0; b < _definition.boundary_conditions.size(); ++b) {
        const auto& bc = _definition.boundary_conditions[b];
        if (bc.type == BoundaryConditionType::dirichlet) {
            std::vector<std::size_t> selected;
            const auto set = std::find_if(source.node_sets().begin(), source.node_sets().end(),
                [&](const NodeSet& s) { return s.name == bc.boundary; });
            if (set != source.node_sets().end())
                selected = set->nodes;
            else {
                for (const auto& side : source.side_set(bc.boundary).sides) {
                    const auto& n = source.elements()[side.element].nodes;
                    selected.insert(
                        selected.end(), {n[side.local_side], n[(side.local_side + 1) % 4], n[4 + side.local_side]});
                }
            }
            std::vector<std::size_t> selected_dofs;
            for (std::size_t r = 0; r < region_count(); ++r)
                for (std::size_t n = 0; n < _meshes[r].nodes().size(); ++n) {
                    if (std::find(selected.begin(), selected.end(), _meshes[r].source_node_ids()[n]) == selected.end())
                        continue;
                    if (bc.field == Field::temperature && !_meshes[r].temperature_nodes()[n]) continue;
                    selected_dofs.push_back(dof(
                        bc.field, bc.field == Field::temperature ? global_temperature_node(r, n) : global_node(r, n)));
                }
            std::sort(selected_dofs.begin(), selected_dofs.end());
            selected_dofs.erase(std::unique(selected_dofs.begin(), selected_dofs.end()), selected_dofs.end());
            if (selected_dofs.empty())
                throw std::invalid_argument("CAX8T Dirichlet boundary contains no active field nodes: " + bc.name);
            for (auto index : selected_dofs) add_dirichlet(index, b);
        } else {
            const auto block = source.side_set_block_id(bc.boundary);
            const auto found = std::find(_block_ids.begin(), _block_ids.end(), block);
            if (found == _block_ids.end())
                throw std::invalid_argument("CAX8T boundary is outside the selected regions");
            const auto r = static_cast<std::size_t>(found - _block_ids.begin());
            for (const auto& edge : _meshes[r].map_side_set(source, bc.boundary).elements)
                _boundaries.push_back({b, r, edge});
            record_configuration_warning(bc, region(r));
        }
    }
    spatial_detail::validate_dirichlet_conditions(
        _dirichlet_conditions, "Conflicting CAX8T Dirichlet conditions", "Duplicate CAX8T Dirichlet conditions");
    refresh_dirichlet_values();
    build_contacts(source);
    _committed_contact_solution = initial_state();
    validate_state(_committed_contact_solution);
}

SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    if (index < volume_contribution_count()) return SpatialContributionType::volume;
    if (index >= volume_contribution_count() + _boundaries.size())
        return _candidates.at(index - volume_contribution_count() - _boundaries.size()).geometry.mechanical
                   ? SpatialContributionType::mechanical_contact
                   : SpatialContributionType::thermal_contact;
    switch (_definition.boundary_conditions[_boundaries.at(index - volume_contribution_count()).definition].type) {
    case BoundaryConditionType::pressure: return SpatialContributionType::pressure;
    case BoundaryConditionType::traction: return SpatialContributionType::traction;
    case BoundaryConditionType::heat_flux: return SpatialContributionType::heat_flux;
    case BoundaryConditionType::convection: return SpatialContributionType::convection;
    default: throw std::logic_error("CAX8T invalid boundary contribution");
    }
}

void SpatialAssembly::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    dofs.clear();
    if (index >= volume_contribution_count() + _boundaries.size()) {
        const auto& values = _candidates.at(index - volume_contribution_count() - _boundaries.size()).dofs;
        dofs.assign(values.begin(), values.end());
        return;
    }
    if (index < volume_contribution_count()) {
        const auto [r, e] = element_location(index);
        const auto& nodes = _meshes[r].elements()[e].nodes;
        for (std::size_t n = 0; n < 4; ++n)
            dofs.push_back(dof(Field::temperature, global_temperature_node(r, nodes[n])));
        for (auto f : {Field::radial_displacement, Field::axial_displacement})
            for (auto n : nodes) dofs.push_back(dof(f, global_node(r, n)));
    } else {
        const auto& b = _boundaries.at(index - volume_contribution_count());
        for (std::size_t n = 0; n < 2; ++n)
            dofs.push_back(dof(Field::temperature, global_temperature_node(b.region, b.edge.nodes[n])));
        for (auto f : {Field::radial_displacement, Field::axial_displacement})
            for (auto n : b.edge.nodes) dofs.push_back(dof(f, global_node(b.region, n)));
    }
}

std::vector<std::size_t> SpatialAssembly::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count()) throw std::out_of_range("CAX8T contribution range");
    std::vector<std::size_t> result, dofs;
    for (auto index = first; index < last; ++index) {
        contribution_dofs(index, dofs);
        result.insert(result.end(), dofs.begin(), dofs.end());
    }
    const auto offset = volume_contribution_count() + _boundaries.size();
    for (const auto& point : _contact_points) {
        if (last <= offset + point.first || first >= offset + point.last) continue;
        for (auto c = point.first; c < point.last; ++c)
            result.insert(result.end(), _candidates[c].dofs.begin(), _candidates[c].dofs.end());
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void SpatialAssembly::validate_local_state(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    if (state.size() != dof_count()) throw std::invalid_argument("CAX8T state layout mismatch");
    const auto dofs = required_state_dofs(first, last);
    for (auto d : dofs)
        if (!std::isfinite(state[d])) throw std::domain_error("CAX8T state must be finite");
    for (auto d : dofs)
        if (d < temperature_node_count() && !(state[d] > 0))
            throw std::domain_error("CAX8T temperature must be positive");
    for (auto index = first; index < std::min(last, volume_contribution_count()); ++index) {
        const auto [r, e] = element_location(index);
        if (region(r).strain_formulation != StrainFormulation::finite) continue;
        auto coordinates = _geometries[r][e].coordinates;
        const auto& nodes = _meshes[r].elements()[e].nodes;
        for (std::size_t n = 0; n < 8; ++n) {
            coordinates[n].r += state[dof(Field::radial_displacement, global_node(r, nodes[n]))];
            coordinates[n].z += state[dof(Field::axial_displacement, global_node(r, nodes[n]))];
        }
        (void)make_quad8_rz_geometry(coordinates);
    }
    validate_contact_state(first, last, state);
}

void SpatialAssembly::compute_boundary(std::size_t index, const std::vector<double>& state,
    std::vector<double>& residual, std::vector<double>* jacobian) const {
    const auto offset = volume_contribution_count() + _boundaries.size();
    if (index >= offset) {
        const auto c = index - offset;
        const auto& candidate = _candidates.at(c);
        if (state.size() != 16) throw std::invalid_argument("CAX8T contact contribution layout mismatch");
        residual.assign(16, 0);
        if (jacobian) jacobian->assign(256, 0);
        if (_active_candidates.at(candidate.point) != c) return;
        std::vector<double> old(16);
        for (std::size_t i = 0; i < 16; ++i) old[i] = _committed_contact_solution[candidate.dofs[i]];
        const auto result =
            compute_line3_contact(candidate.geometry, _heat[candidate.contact], _mechanical[candidate.contact], state,
                old, _contact_histories[candidate.contact][candidate.secondary], jacobian != nullptr);
        residual.assign(result.residual.begin(), result.residual.end());
        if (jacobian) jacobian->assign(result.jacobian.begin(), result.jacobian.end());
        return;
    }
    const auto& boundary = _boundaries.at(index - volume_contribution_count());
    const auto& bc = _definition.boundary_conditions[boundary.definition];
    if (state.size() != 8)
        throw std::invalid_argument("CAX8T quadratic boundary requires eight local degrees of freedom");
    std::array<adlite::Scalar, 8> v;
    for (std::size_t i = 0; i < 8; ++i)
        v[i] = jacobian ? adlite::Scalar::independent(state[i], i, 8) : adlite::Scalar(state[i]);
    std::array<adlite::Scalar, 8> rows{};
    const bool current = boundary_uses_displaced_geometry(bc, region(boundary.region));
    const double load =
        spatial_detail::controlled_value(_definition, _time, _load_factor, bc.value, bc.scale_with_load, bc.function);
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> points = {-g, 0, g}, weights = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    for (std::size_t q = 0; q < 3; ++q) {
        const double x = points[q];
        const std::array<double, 3> shape = {x * (x - 1) / 2, x * (x + 1) / 2, 1 - x * x},
                                    derivative = {x - .5, x + .5, -2 * x};
        const std::array<double, 2> thermal = {(1 - x) / 2, (1 + x) / 2};
        adlite::Scalar radius = 0, tr = 0, tz = 0;
        for (std::size_t n = 0; n < 3; ++n) {
            const auto& point = _meshes[boundary.region].nodes()[boundary.edge.nodes[n]];
            const adlite::Scalar r = current ? point.r + v[2 + n] : adlite::Scalar(point.r),
                                 z = current ? point.z + v[5 + n] : adlite::Scalar(point.z);
            radius += shape[n] * r;
            tr += derivative[n] * r;
            tz += derivative[n] * z;
        }
        const auto length = adlite::hypot(tr, tz), factor = 2 * std::acos(-1.0) * weights[q] * radius,
                   measure = factor * length;
        if (!(length.value() > 0) || !(radius.value() >= 0))
            throw std::domain_error("CAX8T boundary geometry is invalid");
        if (bc.type == BoundaryConditionType::pressure)
            for (std::size_t n = 0; n < 3; ++n) {
                rows[2 + n] += shape[n] * load * factor * tz;
                rows[5 + n] -= shape[n] * load * factor * tr;
            }
        else if (bc.type == BoundaryConditionType::traction) {
            const auto offset = bc.field == Field::radial_displacement ? 2U : 5U;
            for (std::size_t n = 0; n < 3; ++n) rows[offset + n] -= shape[n] * load * measure;
        } else {
            adlite::Scalar flux = -load;
            if (bc.type == BoundaryConditionType::convection) {
                const auto values = convection_values(bc);
                flux = values.coefficient * (thermal[0] * v[0] + thermal[1] * v[1] - values.ambient);
            }
            for (std::size_t n = 0; n < 2; ++n) rows[n] += thermal[n] * flux * measure;
        }
    }
    residual.resize(8);
    if (jacobian) jacobian->assign(64, 0);
    for (std::size_t i = 0; i < 8; ++i) {
        residual[i] = rows[i].value();
        if (jacobian) rows[i].copy_derivatives(jacobian->data() + 8 * i, 8);
    }
}
} // namespace fuelsim::rz8
