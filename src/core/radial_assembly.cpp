#include "radial_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim::radial {
namespace {
constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();
constexpr double pi = 3.141592653589793238462643383279502884;

bool same_axial_coordinate(double a, double b, double height) {
    return std::abs(a - b)
           <= 64.0 * std::numeric_limits<double>::epsilon() * std::max({std::abs(a), std::abs(b), height});
}

void unique_dofs(std::vector<std::size_t>& dofs) {
    std::sort(dofs.begin(), dofs.end());
    dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
}

void validate_history(const ContactPointHistory& history) {
    if (!std::isfinite(history.elastic_tangential_slip) || !std::isfinite(history.total_tangential_slip)
        || history.normal_multiplier != 0.0)
        throw std::invalid_argument("Radial contact history is invalid or uses unsupported contact enforcement");
    for (const auto* values : {&history.cartesian_elastic_tangential_slip,
             &history.cartesian_total_tangential_slip,
             &history.cartesian_contact_normal,
             &history.cartesian_contact_tangent_first})
        for (const auto value : *values)
            if (!std::isfinite(value))
                throw std::invalid_argument("Radial contact history must be finite");
    if (history.cartesian_tangent_basis_initialized
        && (history.cartesian_contact_normal != std::array<double, 3>{{1.0, 0.0, 0.0}}
            || history.cartesian_contact_tangent_first != std::array<double, 3>{{0.0, 0.0, 1.0}}
            || history.cartesian_elastic_tangential_slip
                   != std::array<double, 3>{{0.0, 0.0, history.elastic_tangential_slip}}
            || history.cartesian_total_tangential_slip
                   != std::array<double, 3>{{0.0, 0.0, history.total_tangential_slip}}))
        throw std::invalid_argument("Radial contact history must retain its fixed radial normal and axial tangent");
}
} // namespace

RegionMesh::RegionMesh(const UnstructuredBar2Mesh& source, std::int64_t block)
    : RegionMeshMapping(source, source.nodes().size(), block) {
    std::vector<bool> used(source.nodes().size(), false);
    for (const auto source_element : _source_element_ids)
        for (const auto node : source.elements()[source_element].nodes)
            used[node] = true;
    select_nodes(used);
    for (const auto source_node : _source_node_ids)
        _nodes.push_back(source.nodes()[source_node]);
    for (const auto source_element : _source_element_ids) {
        auto element = source.elements()[source_element];
        for (auto& node : element.nodes)
            node = _source_node_to_local[node];
        _elements.push_back(element);
    }
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredBar2Mesh& source)
    : SpatialLayout(definition,
          spatial_detail::resolve_block_ids(definition, source),
          spatial_detail::DofLayout::axisymmetric_1d),
      _source(source), _source_to_radial(source.nodes().size(), invalid),
      _source_to_axial(source.nodes().size(), invalid), _source_element_to_region(source.elements().size(), invalid) {
    std::vector<std::size_t> node_counts, element_counts, axial_sources;
    std::vector<std::vector<std::size_t>> radial_sources;
    for (std::size_t r = 0; r < region_count(); ++r) {
        if (!region(r).radial_gps)
            throw std::invalid_argument("BAR2 regions require cax2t_gps");
        _meshes.emplace_back(source, _block_ids[r]);
        _materials.emplace_back(region(r).material);
        const auto& mesh = _meshes.back();
        node_counts.push_back(mesh.nodes().size());
        element_counts.push_back(mesh.elements().size());
        radial_sources.push_back(mesh.source_node_ids());
        _geometries.emplace_back();
        for (std::size_t e = 0; e < mesh.elements().size(); ++e) {
            const auto& element = mesh.elements()[e];
            _source_element_to_region[mesh.source_element_ids()[e]] = r;
            _geometries.back().push_back(elements::make_cax2t_gps_geometry(
                {{mesh.nodes()[element.nodes[0]].r, mesh.nodes()[element.nodes[1]].r}},
                source.nodes()[element.axial_nodes[0]].z,
                source.nodes()[element.axial_nodes[1]].z));
            axial_sources.insert(axial_sources.end(), element.axial_nodes.begin(), element.axial_nodes.end());
        }
    }
    unique_dofs(axial_sources);
    initialize_counts(node_counts, element_counts);
    initialize_radial_shared_nodes(radial_sources, axial_sources.size());
    for (std::size_t a = 0; a < axial_sources.size(); ++a)
        _source_to_axial[axial_sources[a]] = a;
    for (std::size_t r = 0; r < region_count(); ++r)
        for (std::size_t n = 0; n < radial_sources[r].size(); ++n)
            _source_to_radial[radial_sources[r][n]] = global_node(r, n);
    refresh_heat_sources();
    build_boundaries();
    build_contacts();
    _committed_contact_solution = initial_state();
    validate_state(_committed_contact_solution);
}

std::size_t SpatialAssembly::axial_dof(std::size_t source_node) const {
    if (source_node >= _source_to_axial.size() || _source_to_axial[source_node] == invalid)
        throw std::invalid_argument("Radial source node does not carry an active axial control field");
    return dof(Field::axial_displacement, _source_to_axial[source_node]);
}

std::size_t SpatialAssembly::radial_node(std::size_t source_node) const {
    if (source_node >= _source_to_radial.size() || _source_to_radial[source_node] == invalid)
        throw std::invalid_argument("Radial source node does not carry active temperature and radial fields");
    return _source_to_radial[source_node];
}

std::vector<double> SpatialAssembly::reference_state() const {
    std::vector<double> state(dof_count(), 0.0);
    for (std::size_t r = 0; r < region_count(); ++r)
        for (std::size_t n = 0; n < region_mesh(r).nodes().size(); ++n)
            state[dof(Field::temperature, global_temperature_node(r, n))] = region(r).initial_temperature;
    return state;
}

std::array<std::size_t, 2> SpatialAssembly::region_axial_dofs(std::size_t r, std::size_t e) const {
    const auto& nodes = region_mesh(r).elements().at(e).axial_nodes;
    return {{axial_dof(nodes[0]), axial_dof(nodes[1])}};
}

void SpatialAssembly::refresh_heat_sources() {
    _heat_sources.resize(region_count());
    for (std::size_t r = 0; r < region_count(); ++r)
        _heat_sources[r] = region_heat_source(r);
}

void SpatialAssembly::set_load_factor(double value) {
    set_load_factor_value(value);
    refresh_dirichlet_values();
    refresh_heat_sources();
}

void SpatialAssembly::set_time(double value) {
    set_time_value(value);
    refresh_dirichlet_values();
    refresh_heat_sources();
}

void SpatialAssembly::set_heat_source_interval(double begin, double end) {
    if (!std::isfinite(begin) || !std::isfinite(end) || !(end > begin))
        throw std::invalid_argument("Radial heat-source interval must be finite and increasing");
    for (std::size_t r = 0; r < region_count(); ++r)
        _heat_sources[r] = region_heat_source_average(r, begin, end);
}

std::vector<ElementSide> SpatialAssembly::boundary_sides(const std::string& name) const {
    const auto& sides = _source.side_set(name).sides;
    if (sides.empty())
        throw std::invalid_argument("Radial boundary side set is empty: " + name);
    std::vector<std::pair<std::size_t, std::size_t>> seen;
    for (const auto& side : sides) {
        if (_source_element_to_region.at(side.element) == invalid)
            throw std::invalid_argument("Radial boundary belongs to an unselected region: " + name);
        const auto key = std::make_pair(side.element, side.local_side);
        if (std::find(seen.begin(), seen.end(), key) != seen.end())
            throw std::invalid_argument("Radial boundary contains a repeated endpoint: " + name);
        seen.push_back(key);
    }
    return sides;
}

void SpatialAssembly::build_boundaries() {
    for (std::size_t b = 0; b < _definition.boundary_conditions.size(); ++b) {
        const auto& bc = _definition.boundary_conditions[b];
        if (bc.name.empty() || bc.boundary.empty() || !std::isfinite(bc.value)
            || (bc.scale_with_load && !bc.function.empty()))
            throw std::invalid_argument("Radial boundary names, values, and controls are invalid");
        const auto set = std::find_if(_source.node_sets().begin(), _source.node_sets().end(), [&](const NodeSet& item) {
            return item.name == bc.boundary;
        });
        const bool axial = bc.type == BoundaryConditionType::axial_force
                           || (bc.type == BoundaryConditionType::dirichlet && bc.field == Field::axial_displacement);
        if (axial) {
            if (set == _source.node_sets().end() || set->nodes.empty())
                throw std::invalid_argument("Radial axial constraints and forces require an axial control node set");
            std::vector<std::size_t> selected = set->nodes;
            unique_dofs(selected);
            if (selected.size() != set->nodes.size())
                throw std::invalid_argument("Radial axial node set contains duplicate nodes");
            for (const auto source : selected) {
                const auto index = axial_dof(source);
                if (bc.type == BoundaryConditionType::dirichlet)
                    add_dirichlet(index, b);
                else
                    _boundaries.push_back({b, invalid, invalid, source});
            }
            continue;
        }
        if (bc.type == BoundaryConditionType::dirichlet) {
            if (bc.field != Field::temperature && bc.field != Field::radial_displacement)
                throw std::invalid_argument(
                    "Radial Dirichlet field must be temperature, radial, or axial displacement");
            std::vector<std::size_t> selected;
            if (set != _source.node_sets().end())
                selected = set->nodes;
            else
                for (const auto& side : boundary_sides(bc.boundary))
                    selected.push_back(_source.elements()[side.element].nodes[side.local_side]);
            unique_dofs(selected);
            if (selected.empty())
                throw std::invalid_argument("Radial Dirichlet boundary contains no nodes");
            for (const auto source : selected) {
                if (bc.field == Field::radial_displacement && _source.nodes()[source].r == 0.0 && bc.value != 0.0)
                    throw std::invalid_argument("Radial displacement at the symmetry axis must be zero");
                add_dirichlet(dof(bc.field, radial_node(source)), b);
            }
            continue;
        }
        if (bc.type != BoundaryConditionType::pressure && bc.type != BoundaryConditionType::traction
            && bc.type != BoundaryConditionType::heat_flux && bc.type != BoundaryConditionType::convection)
            throw std::invalid_argument("Radial boundary type is unsupported");
        if (bc.type == BoundaryConditionType::traction && bc.field != Field::radial_displacement
            && bc.field != Field::axial_displacement)
            throw std::invalid_argument("Radial cylindrical traction requires radial or axial displacement field");
        if (bc.type == BoundaryConditionType::convection
            && (!std::isfinite(bc.heat_transfer_coefficient) || bc.heat_transfer_coefficient < 0.0
                || !std::isfinite(bc.ambient_temperature) || !(bc.ambient_temperature > 0.0)))
            throw std::invalid_argument("Radial convection requires nonnegative coefficient and positive temperature");
        for (const auto& side : boundary_sides(bc.boundary)) {
            _boundaries.push_back({b, side.element, side.local_side, invalid});
            record_configuration_warning(bc, region(_source_element_to_region[side.element]));
        }
    }
    for (const auto source : _source.radial_source_node_ids()) {
        if (_source_to_radial[source] == invalid || _source.nodes()[source].r != 0.0)
            continue;
        const auto axis = dof(Field::radial_displacement, radial_node(source));
        const auto condition = std::find_if(_dirichlet_conditions.begin(),
            _dirichlet_conditions.end(),
            [&](const DirichletCondition& item) { return item.dof == axis; });
        if (condition == _dirichlet_conditions.end())
            _dirichlet_conditions.push_back({axis, 0.0});
        else if (condition->value != 0.0)
            throw std::invalid_argument("Radial displacement at the symmetry axis must be zero");
    }
    spatial_detail::validate_dirichlet_conditions(_dirichlet_conditions,
        "Conflicting radial Dirichlet conditions",
        "Duplicate radial Dirichlet conditions");
    refresh_dirichlet_values();
}

void SpatialAssembly::build_contacts() {
    _contact_histories.resize(_definition.contacts.size());
    std::vector<bool> thermal_secondary(_source.nodes().size(), false),
        mechanical_secondary(_source.nodes().size(), false);
    for (std::size_t c = 0; c < _definition.contacts.size(); ++c) {
        const auto& definition = _definition.contacts[c];
        if (definition.mechanical_formulation != MechanicalContactFormulation::penalty || definition.automatic_penalty
            || definition.mechanical_sliding != MechanicalContactSliding::small
            || definition.mechanical_discretization == MechanicalContactDiscretization::node_to_surface
            || definition.thermal_discretization != ThermalContactDiscretization::surface_to_surface)
            throw std::invalid_argument(
                "Radial contact requires explicit penalty, small sliding, and surface-to-surface discretization");
        GapHeatProperties heat{definition.gap_conductivity, definition.minimum_gap};
        heat.law = definition.gap_heat_conductance_law;
        heat.conductance = definition.gap_conductance;
        heat.clearance_derivative = definition.gap_conductance_clearance_derivative;
        heat.pressure_derivative = definition.gap_conductance_pressure_derivative;
        heat.temperature_derivative = definition.gap_conductance_temperature_derivative;
        heat.reference_temperature = definition.gap_conductance_reference_temperature;
        heat.contact_penalty = definition.mechanical ? definition.penalty : 0.0;
        _heat.push_back(heat);
        const auto primary = boundary_sides(definition.primary), secondary = boundary_sides(definition.secondary);
        if (primary.size() != secondary.size())
            throw std::invalid_argument("Radial contact requires matching one-to-one axial segmentation");
        double average_primary_length = 0.0;
        for (const auto& side : primary) {
            const auto& element = _source.elements()[side.element];
            average_primary_length +=
                (_source.nodes()[element.axial_nodes[1]].z - _source.nodes()[element.axial_nodes[0]].z)
                / static_cast<double>(primary.size());
        }
        const double elastic_slip = definition.friction_slip_tolerance * average_primary_length;
        if (!std::isfinite(elastic_slip))
            throw std::invalid_argument("Radial contact elastic slip length must be finite");
        _normal.push_back(
            {definition.mechanical ? definition.penalty : 0.0, definition.friction_coefficient, false, elastic_slip});
        std::vector<bool> used(primary.size(), false);
        for (const auto& s : secondary) {
            if (s.local_side != 1)
                throw std::invalid_argument("Radial contact secondary must be the outer side of the inner cylinder");
            const auto& se = _source.elements()[s.element];
            const double lower = _source.nodes()[se.axial_nodes[0]].z, upper = _source.nodes()[se.axial_nodes[1]].z;
            std::size_t selected = invalid;
            for (std::size_t p = 0; p < primary.size(); ++p) {
                if (primary[p].local_side != 0)
                    throw std::invalid_argument("Radial contact primary must be the inner side of the outer cylinder");
                const auto& pe = _source.elements()[primary[p].element];
                if (same_axial_coordinate(lower, _source.nodes()[pe.axial_nodes[0]].z, upper - lower)
                    && same_axial_coordinate(upper, _source.nodes()[pe.axial_nodes[1]].z, upper - lower)) {
                    if (selected != invalid)
                        throw std::invalid_argument("Radial contact axial segment has multiple primary candidates");
                    selected = p;
                }
            }
            if (selected == invalid || used[selected])
                throw std::invalid_argument("Radial contact requires matching unique axial segments");
            used[selected] = true;
            const auto& pe = _source.elements()[primary[selected].element];
            const auto formulation = region(_source_element_to_region[s.element]).strain_formulation;
            if (formulation != region(_source_element_to_region[primary[selected].element]).strain_formulation)
                throw std::invalid_argument("Radial contact cannot mix small and finite strain regions");
            const auto pn = pe.nodes[0], sn = se.nodes[1];
            if (pn == sn || !(_source.nodes()[pn].r > 0.0) || !(_source.nodes()[sn].r > 0.0))
                throw std::invalid_argument("Radial contact requires distinct positive-radius surface nodes");
            if ((definition.thermal && thermal_secondary[sn]) || (definition.mechanical && mechanical_secondary[sn]))
                throw std::invalid_argument(
                    "Radial secondary source endpoint has duplicate contact constraints through boundary aliases");
            thermal_secondary[sn] = thermal_secondary[sn] || definition.thermal;
            mechanical_secondary[sn] = mechanical_secondary[sn] || definition.mechanical;
            Contact contact{{_source.nodes()[pn].r, _source.nodes()[sn].r, lower, upper},
                {{dof(Field::temperature, radial_node(pn)),
                    dof(Field::temperature, radial_node(sn)),
                    dof(Field::radial_displacement, radial_node(pn)),
                    dof(Field::radial_displacement, radial_node(sn)),
                    axial_dof(pe.axial_nodes[0]),
                    axial_dof(pe.axial_nodes[1]),
                    axial_dof(se.axial_nodes[0]),
                    axial_dof(se.axial_nodes[1])}},
                c,
                _contact_histories[c].size(),
                pn,
                sn,
                primary[selected].element,
                s.element,
                formulation};
            _contacts.push_back(contact);
            _contact_histories[c].resize(_contact_histories[c].size() + ring_gps_quadrature_point_count);
        }
    }
}

std::pair<std::size_t, std::size_t> SpatialAssembly::contribution_partition(std::size_t rank, std::size_t ranks) const {
    if (ranks == 0 || rank >= ranks)
        throw std::invalid_argument("Radial contribution partition rank is invalid");
    const auto count = contribution_count(), chunk = count / ranks, remainder = count % ranks;
    const auto begin = rank * chunk + std::min(rank, remainder);
    return {begin, begin + chunk + (rank < remainder ? 1U : 0U)};
}

SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    if (index < volume_contribution_count())
        return SpatialContributionType::volume;
    index -= volume_contribution_count();
    if (index < _contacts.size())
        return _definition.contacts[_contacts[index].definition].mechanical
                   ? SpatialContributionType::mechanical_contact
                   : SpatialContributionType::thermal_contact;
    const auto& boundary = _boundaries.at(index - _contacts.size());
    switch (_definition.boundary_conditions[boundary.definition].type) {
    case BoundaryConditionType::pressure:
        return SpatialContributionType::pressure;
    case BoundaryConditionType::traction:
    case BoundaryConditionType::axial_force:
        return SpatialContributionType::traction;
    case BoundaryConditionType::heat_flux:
        return SpatialContributionType::heat_flux;
    case BoundaryConditionType::convection:
        return SpatialContributionType::convection;
    default:
        throw std::logic_error("Radial boundary contribution type is invalid");
    }
}

void SpatialAssembly::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    dofs.clear();
    if (index < volume_contribution_count()) {
        const auto [r, e] = element_location(index);
        const auto& element = region_mesh(r).elements()[e];
        for (const auto field : {Field::temperature, Field::radial_displacement})
            for (const auto n : element.nodes)
                dofs.push_back(dof(field, global_node(r, n)));
        for (const auto source : element.axial_nodes)
            dofs.push_back(axial_dof(source));
        return;
    }
    index -= volume_contribution_count();
    if (index < _contacts.size()) {
        const auto& contact = _contacts[index];
        dofs.assign(contact.dofs.begin(), contact.dofs.end());
        unique_dofs(dofs);
        return;
    }
    const auto& boundary = _boundaries.at(index - _contacts.size());
    if (boundary.axial_source != invalid) {
        dofs.push_back(axial_dof(boundary.axial_source));
        return;
    }
    const auto& element = _source.elements()[boundary.element];
    const auto node = radial_node(element.nodes[boundary.side]);
    dofs = {dof(Field::temperature, node),
        dof(Field::radial_displacement, node),
        axial_dof(element.axial_nodes[0]),
        axial_dof(element.axial_nodes[1])};
}

Cax2tGpsLocalValues SpatialAssembly::volume_state(std::size_t index, const std::vector<double>& global) const {
    if (index >= volume_contribution_count() || global.size() != dof_count())
        throw std::invalid_argument("Radial volume state layout or index is invalid");
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    Cax2tGpsLocalValues values{};
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] = global[dofs[i]];
    return values;
}

std::vector<std::size_t> SpatialAssembly::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("Radial contribution range is invalid");
    std::vector<std::size_t> result, dofs;
    for (auto i = first; i < last; ++i) {
        contribution_dofs(i, dofs);
        result.insert(result.end(), dofs.begin(), dofs.end());
    }
    unique_dofs(result);
    return result;
}

void SpatialAssembly::validate_local_state(std::size_t first,
    std::size_t last,
    const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Radial global state size is invalid");
    for (const auto d : required_state_dofs(first, last))
        if (!std::isfinite(state[d]) || (d < temperature_node_count() && !(state[d] > 0.0)))
            throw std::domain_error("Radial state must be finite with positive temperatures");
    for (auto i = first; i < std::min(last, volume_contribution_count()); ++i) {
        const auto [r, e] = element_location(i);
        const auto local = volume_state(i, state);
        const auto& geometry = region_element_geometry(r, e);
        const double inner = geometry.radii[0] + local[2], outer = geometry.radii[1] + local[3];
        const double height = geometry.z_upper - geometry.z_lower + local[5] - local[4];
        if ((geometry.radii[0] == 0.0 && local[2] != 0.0) || !std::isfinite(inner) || !std::isfinite(outer)
            || !std::isfinite(height) || inner < 0.0 || !(outer > inner) || !(height > 0.0))
            throw std::domain_error(
                "Radial trial geometry must preserve the axis, ordered radii, and positive axial height");
    }
    const auto offset = volume_contribution_count();
    for (auto i = std::max(first, offset); i < std::min(last, offset + _contacts.size()); ++i)
        (void)evaluate_global_contact(i - offset, state, false);
}

elements::Cax2tGpsResult SpatialAssembly::evaluate_volume(std::size_t r,
    std::size_t e,
    const std::vector<double>& global,
    const std::vector<double>& committed,
    const Cax2tGpsMaterialHistory* history,
    double dt,
    double time,
    bool include_thermal,
    bool linearize) const {
    const auto index = region_element_offset(r) + e;
    const auto state = volume_state(index, global);
    const auto old =
        committed.empty()
            ? Cax2tGpsLocalValues{{region(r).initial_temperature, region(r).initial_temperature, 0.0, 0.0, 0.0, 0.0}}
            : volume_state(index, committed);
    const elements::Cax2tGpsInput input{_materials.at(r),
        region_element_geometry(r, e),
        state,
        old,
        history,
        dt,
        time,
        _heat_sources.at(r),
        region(r).strain_formulation,
        include_thermal};
    return elements::evaluate_cax2t_gps(input, {true, linearize, true, true});
}

elements::RingGpsResult SpatialAssembly::evaluate_contact(std::size_t index,
    const RingGpsLocalValues& state,
    bool linearize,
    bool history) const {
    const auto& contact = _contacts.at(index);
    const auto& definition = _definition.contacts[contact.definition];
    RingGpsLocalValues old{};
    for (std::size_t i = 0; i < old.size(); ++i)
        old[i] = _committed_contact_solution.at(contact.dofs[i]);
    RingGpsHistory committed{};
    for (std::size_t q = 0; q < committed.size(); ++q)
        committed[q] = _contact_histories[contact.definition][contact.history_offset + q];
    return elements::evaluate_ring_gps({contact.geometry,
                                           state,
                                           old,
                                           committed,
                                           definition.thermal,
                                           definition.mechanical,
                                           _heat[contact.definition],
                                           _normal[contact.definition],
                                           contact.strain_formulation},
        {true, linearize, history, false});
}

elements::RingGpsResult
SpatialAssembly::evaluate_global_contact(std::size_t index, const std::vector<double>& global, bool history) const {
    if (global.size() != dof_count())
        throw std::invalid_argument("Radial contact global state size is invalid");
    const auto& contact = _contacts.at(index);
    RingGpsLocalValues state{};
    for (std::size_t i = 0; i < state.size(); ++i)
        state[i] = global[contact.dofs[i]];
    return evaluate_contact(index, state, false, history);
}

void SpatialAssembly::compute_contribution(std::size_t index,
    const std::vector<double>& local_state,
    const std::vector<double>& committed_global,
    const Cax2tGpsMaterialHistory* history,
    double dt,
    double time,
    bool include_thermal,
    bool linearize,
    LocalContribution& result) const {
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    if (local_state.size() != dofs.size())
        throw std::invalid_argument("Radial local contribution state size is invalid");
    result = {};
    if (index < volume_contribution_count()) {
        const auto [r, e] = element_location(index);
        Cax2tGpsLocalValues state{};
        std::copy(local_state.begin(), local_state.end(), state.begin());
        const auto old = committed_global.empty() ? Cax2tGpsLocalValues{{region(r).initial_temperature,
                                                        region(r).initial_temperature,
                                                        0.0,
                                                        0.0,
                                                        0.0,
                                                        0.0}}
                                                  : volume_state(index, committed_global);
        result.volume = elements::evaluate_cax2t_gps({_materials[r],
                                                         _geometries[r][e],
                                                         state,
                                                         old,
                                                         history,
                                                         dt,
                                                         time,
                                                         _heat_sources[r],
                                                         region(r).strain_formulation,
                                                         include_thermal},
            {true, linearize, true, true});
        result.residual.assign(result.volume.residual.begin(), result.volume.residual.end());
        if (linearize)
            result.jacobian.assign(result.volume.jacobian.begin(), result.volume.jacobian.end());
    } else if (index < volume_contribution_count() + _contacts.size()) {
        const auto c = index - volume_contribution_count();
        const auto& contact = _contacts[c];
        RingGpsLocalValues state{};
        std::array<std::size_t, ring_gps_local_dof_count> mapping{};
        for (std::size_t i = 0; i < state.size(); ++i) {
            mapping[i] =
                static_cast<std::size_t>(std::lower_bound(dofs.begin(), dofs.end(), contact.dofs[i]) - dofs.begin());
            state[i] = local_state[mapping[i]];
        }
        result.contact = evaluate_contact(c, state, linearize);
        result.residual.assign(dofs.size(), 0.0);
        if (linearize)
            result.jacobian.assign(dofs.size() * dofs.size(), 0.0);
        for (std::size_t i = 0; i < state.size(); ++i) {
            result.residual[mapping[i]] += result.contact.residual[i];
            if (linearize)
                for (std::size_t j = 0; j < state.size(); ++j)
                    result.jacobian[mapping[i] * dofs.size() + mapping[j]] +=
                        result.contact.jacobian[i * state.size() + j];
        }
    } else
        compute_surface(index - volume_contribution_count() - _contacts.size(), local_state, linearize, result);
}

void SpatialAssembly::compute_surface(std::size_t index,
    const std::vector<double>& state,
    bool linearize,
    LocalContribution& result) const {
    const auto& boundary = _boundaries.at(index);
    const auto& bc = _definition.boundary_conditions[boundary.definition];
    const double load =
        spatial_detail::controlled_value(_definition, _time, _load_factor, bc.value, bc.scale_with_load, bc.function);
    if (!std::isfinite(load))
        throw std::domain_error("Radial boundary value must be finite");
    result.residual.assign(state.size(), 0.0);
    if (linearize)
        result.jacobian.assign(state.size() * state.size(), 0.0);
    if (boundary.axial_source != invalid) {
        result.residual[0] = -load;
        return;
    }
    const auto& element = _source.elements()[boundary.element];
    const bool displaced = boundary_uses_displaced_geometry(bc, region(_source_element_to_region[boundary.element]));
    const double radius = _source.nodes()[element.nodes[boundary.side]].r + (displaced ? state[1] : 0.0);
    const double height = _source.nodes()[element.axial_nodes[1]].z - _source.nodes()[element.axial_nodes[0]].z
                          + (displaced ? state[3] - state[2] : 0.0);
    const double area = 2.0 * pi * radius * height;
    if (!std::isfinite(area) || radius < 0.0 || !(height > 0.0))
        throw std::domain_error("Radial cylindrical boundary requires finite nonnegative radius and positive height");
    const std::array<double, 4> area_derivative =
        displaced ? std::array<double, 4>{{0.0, 2.0 * pi * height, -2.0 * pi * radius, 2.0 * pi * radius}}
                  : std::array<double, 4>{};
    std::array<double, 4> density{};
    switch (bc.type) {
    case BoundaryConditionType::pressure:
        density[1] = load * (boundary.side == 0 ? -1.0 : 1.0);
        break;
    case BoundaryConditionType::traction:
        if (bc.field == Field::radial_displacement)
            density[1] = -load;
        else {
            density[2] = -0.5 * load;
            density[3] = -0.5 * load;
        }
        break;
    case BoundaryConditionType::heat_flux:
        result.surface_heat_input_rate = load * area;
        density[0] = -load;
        break;
    case BoundaryConditionType::convection: {
        const auto values = convection_values(bc);
        if (!std::isfinite(values.coefficient) || values.coefficient < 0.0 || !std::isfinite(values.ambient)
            || !(values.ambient > 0.0))
            throw std::domain_error("Radial convection requires nonnegative coefficient and positive temperature");
        result.convection_heat_rate = values.coefficient * area * (state[0] - values.ambient);
        density[0] = values.coefficient * (state[0] - values.ambient);
        if (linearize)
            result.jacobian[0] = values.coefficient * area;
        break;
    }
    default:
        throw std::logic_error("Radial surface boundary type is invalid");
    }
    for (std::size_t i = 0; i < density.size(); ++i) {
        result.residual[i] = density[i] * area;
        if (linearize)
            for (std::size_t j = 0; j < area_derivative.size(); ++j)
                result.jacobian[i * density.size() + j] += density[i] * area_derivative[j];
    }
}

void SpatialAssembly::compute_boundary(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian) const {
    if (index < volume_contribution_count())
        throw std::invalid_argument("Radial compute_boundary cannot evaluate a volume contribution");
    LocalContribution result;
    compute_contribution(index,
        state,
        _committed_contact_solution,
        nullptr,
        0.0,
        _time,
        false,
        jacobian != nullptr,
        result);
    residual = std::move(result.residual);
    if (jacobian)
        *jacobian = std::move(result.jacobian);
}

double SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    validate_state(state);
    auto histories = _contact_histories;
    std::vector<double> solution = state;
    double dissipation = 0.0;
    for (std::size_t c = 0; c < _contacts.size(); ++c) {
        const auto result = evaluate_global_contact(c, state);
        const auto& contact = _contacts[c];
        for (std::size_t q = 0; q < result.history.size(); ++q)
            histories[contact.definition][contact.history_offset + q] = result.history[q];
        dissipation += result.friction_dissipation;
    }
    if (!std::isfinite(dissipation))
        throw std::domain_error("Radial contact friction dissipation must be finite");
    _contact_histories.swap(histories);
    _committed_contact_solution.swap(solution);
    return dissipation;
}

void SpatialAssembly::restore_contact_state(const std::vector<double>& state,
    std::vector<std::vector<ContactPointHistory>> histories) {
    if (state.size() != dof_count() || histories.size() != _contact_histories.size())
        throw std::invalid_argument("Radial restored contact state layout is invalid");
    for (std::size_t c = 0; c < histories.size(); ++c) {
        if (histories[c].size() != _contact_histories[c].size())
            throw std::invalid_argument("Radial restored contact history count is invalid");
        for (const auto& history : histories[c])
            validate_history(history);
    }
    validate_state(state);
    std::vector<double> solution = state;
    _contact_histories.swap(histories);
    _committed_contact_solution.swap(solution);
}

InterfaceSummary SpatialAssembly::summarize_interface(std::size_t contact_id, const std::vector<double>& state) const {
    (void)_definition.contacts.at(contact_id);
    InterfaceSummary summary;
    for (std::size_t c = 0; c < _contacts.size(); ++c) {
        if (_contacts[c].definition != contact_id)
            continue;
        const auto result = evaluate_global_contact(c, state, false);
        summary.total_heat_rate += result.heat_rate;
        summary.total_contact_force += result.normal_force;
        summary.total_tangential_force += result.tangential_force;
        for (const auto& point : result.points) {
            summary.minimum_gap = std::min(summary.minimum_gap, point.gap);
            summary.minimum_contact_gap = std::min(summary.minimum_contact_gap, point.gap);
            summary.maximum_contact_pressure = std::max(summary.maximum_contact_pressure, point.pressure);
            ++summary.projected_contact_nodes;
            if (point.pressure > 0.0)
                ++summary.active_contact_nodes;
        }
    }
    return summary;
}

std::vector<ContactNodeSummary> SpatialAssembly::summarize_contact_nodes(std::size_t contact_id,
    const std::vector<double>& state) const {
    (void)_definition.contacts.at(contact_id);
    std::vector<ContactNodeSummary> rows;
    constexpr double gauss = 0.57735026918962576451;
    for (std::size_t c = 0; c < _contacts.size(); ++c) {
        const auto& contact = _contacts[c];
        if (contact.definition != contact_id)
            continue;
        const auto result = evaluate_global_contact(c, state, false);
        const double height = contact.geometry.z_upper - contact.geometry.z_lower;
        for (std::size_t q = 0; q < result.points.size(); ++q) {
            const auto& point = result.points[q];
            const double shape = 0.5 * (1.0 + (q == 0 ? -gauss : gauss));
            rows.push_back({contact.geometry.secondary_reference_radius,
                contact.geometry.z_lower + shape * height,
                true,
                contact.primary_source,
                point.gap,
                point.pressure,
                point.weighted_measure,
                0.5 * height,
                point.pressure * point.weighted_measure,
                point.signed_tangential_traction,
                point.signed_tangential_traction * point.weighted_measure,
                point.elastic_tangential_slip,
                point.sliding,
                point.total_tangential_slip});
        }
    }
    return rows;
}

std::vector<std::size_t> SpatialAssembly::contact_secondary_source_nodes(std::size_t contact_id) const {
    (void)_definition.contacts.at(contact_id);
    std::vector<std::size_t> nodes;
    for (const auto& contact : _contacts)
        if (contact.definition == contact_id)
            for (std::size_t q = 0; q < ring_gps_quadrature_point_count; ++q)
                nodes.push_back(contact.secondary_source);
    return nodes;
}

std::vector<std::pair<std::size_t, std::size_t>> SpatialAssembly::contact_source_elements(
    std::size_t contact_id) const {
    (void)_definition.contacts.at(contact_id);
    std::vector<std::pair<std::size_t, std::size_t>> result;
    for (const auto& contact : _contacts)
        if (contact.definition == contact_id)
            result.emplace_back(contact.primary_element, contact.secondary_element);
    return result;
}

elements::RingGpsResult
SpatialAssembly::evaluate_contact(std::size_t contact_id, std::size_t pair, const std::vector<double>& state) const {
    (void)_definition.contacts.at(contact_id);
    std::size_t selected = 0;
    for (std::size_t i = 0; i < _contacts.size(); ++i)
        if (_contacts[i].definition == contact_id) {
            if (selected == pair)
                return evaluate_global_contact(i, state);
            ++selected;
        }
    throw std::out_of_range("Radial contact pair index is invalid");
}
} // namespace fuelsim::radial
