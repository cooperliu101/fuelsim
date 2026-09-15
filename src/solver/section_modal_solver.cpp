#include "solver/section_modal_solver.hpp"
#include "modal_condensation.hpp"
#include "quad8_face.hpp"
#include "solver/reduced_section_basis.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <petscksp.h>
#include <set>
#include <stdexcept>
#include <string>

namespace fuelsim {
namespace {
void check(PetscErrorCode code) {
    if (code != PETSC_SUCCESS)
        throw std::runtime_error("Modal axial PETSc failure: " + std::to_string(code));
}

std::vector<double> gather(Vec vector) {
    Vec all = nullptr;
    VecScatter scatter = nullptr;
    check(VecScatterCreateToAll(vector, &scatter, &all));
    check(VecScatterBegin(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD));
    check(VecScatterEnd(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD));
    PetscInt n;
    check(VecGetSize(all, &n));
    const PetscScalar* values;
    check(VecGetArrayRead(all, &values));
    std::vector<double> result(values, values + n);
    check(VecRestoreArrayRead(all, &values));
    check(VecScatterDestroy(&scatter));
    check(VecDestroy(&all));
    return result;
}

struct Objects final {
    Mat stiffness = nullptr, constrained = nullptr;
    Vec diagonal = nullptr, physical_load = nullptr, rhs = nullptr, solution = nullptr, defect = nullptr;
    Vec correction = nullptr;
    KSP solver = nullptr;

    ~Objects() {
        (void)KSPDestroy(&solver);
        (void)VecDestroy(&diagonal);
        (void)VecDestroy(&physical_load);
        (void)VecDestroy(&rhs);
        (void)VecDestroy(&solution);
        (void)VecDestroy(&defect);
        (void)VecDestroy(&correction);
        (void)MatDestroy(&stiffness);
        (void)MatDestroy(&constrained);
    }
};

struct AxialLayout final {
    std::vector<std::size_t> modes;
    std::vector<PetscInt> indices;
    std::vector<bool> eliminated;
    std::size_t count = 0, retained = 0;
};

struct Geometry final {
    std::unique_ptr<CrossSection> section;
    std::vector<double> axial;
    std::vector<std::size_t> source_to_section, element_section;
    AxialLayout layout;
};

Geometry extract(const UnstructuredHex20Mesh& source, const SpatialDefinition& definition) {
    if (!definition.contacts.empty() || !definition.time_tables.empty())
        throw std::invalid_argument("Section modal mechanics does not support contact or time functions");
    Geometry result;
    std::set<double> axial;
    std::vector<SectionRegion> regions;
    std::map<std::int64_t, std::size_t> block_regions;
    for (const auto& region : definition.regions) {
        if (region.requested_cartesian_node_count != 20
            || region.hex20_element_formulation != Hex20ElementFormulation::c3d20t
            || region.strain_formulation != StrainFormulation::small || region.volumetric_heat_source != 0.0
            || region.body_acceleration != std::array<double, 3>{})
            throw std::invalid_argument("Section modal regions require small strain and no thermal or body sources");
        const auto id = source.element_block(region.block).id;
        if (!block_regions.emplace(id, regions.size()).second)
            throw std::invalid_argument("Duplicate section modal material block");
        regions.push_back({region.name, IsotropicThermoelasticMaterial(region.material), region.initial_temperature});
    }
    const double z0 = std::min_element(source.nodes().begin(), source.nodes().end(), [](const auto& a, const auto& b) {
        return a.z < b.z;
    })->z;
    std::map<std::pair<double, double>, std::size_t> xy;
    std::vector<CartesianPoint3> nodes;
    for (const auto& point : source.nodes())
        if (point.z == z0) {
            if (!xy.emplace(std::make_pair(point.x, point.y), nodes.size()).second)
                throw std::invalid_argument("Section modal mesh has duplicate unconnected cross-section coordinates");
            nodes.push_back({point.x, point.y, 0.0});
        }
    for (const auto& node : source.nodes()) {
        const auto found = xy.find({node.x, node.y});
        if (found == xy.end())
            throw std::invalid_argument("Section modal mesh is not a fixed extrusion along z");
        result.source_to_section.push_back(found->second);
    }
    constexpr std::array<std::size_t, 8> bottom{0, 1, 2, 3, 8, 9, 10, 11}, top{4, 5, 6, 7, 16, 17, 18, 19};
    std::vector<SectionCell> cells;
    std::map<std::array<std::size_t, 9>, std::size_t> sections;
    for (const auto& element : source.elements()) {
        const double lower = source.nodes()[element.nodes[0]].z, upper = source.nodes()[element.nodes[4]].z;
        if (!(upper > lower))
            throw std::invalid_argument("Section modal HEX20 orientation must follow positive z");
        axial.insert(lower);
        axial.insert(upper);
        for (std::size_t i = 0; i < 8; ++i) {
            const auto lo = element.nodes[bottom[i]], hi = element.nodes[top[i]];
            if (source.nodes()[lo].z != lower || source.nodes()[hi].z != upper
                || result.source_to_section[lo] != result.source_to_section[hi])
                throw std::invalid_argument("Section modal HEX20 end faces do not match");
        }
        for (std::size_t i = 0; i < 4; ++i)
            if (result.source_to_section[element.nodes[12 + i]] != result.source_to_section[element.nodes[i]]
                || std::abs(source.nodes()[element.nodes[12 + i]].z - 0.5 * (lower + upper))
                       > 1.0e-12 * (upper - lower))
                throw std::invalid_argument("Section modal axial edge geometry must be straight");
    }
    result.axial.assign(axial.begin(), axial.end());
    for (std::size_t e = 0; e < source.elements().size(); ++e) {
        if (source.nodes()[source.elements()[e].nodes[0]].z != z0)
            continue;
        const auto found = block_regions.find(source.element_block_ids()[e]);
        if (found == block_regions.end())
            throw std::invalid_argument("Section modal mesh contains an unselected block");
        SectionCell cell{};
        cell.region = found->second;
        std::array<std::size_t, 9> key{};
        for (std::size_t i = 0; i < 8; ++i)
            key[i] = cell.element.nodes[i] = result.source_to_section[source.elements()[e].nodes[bottom[i]]];
        key[8] = cell.region;
        if (!sections.emplace(key, cells.size()).second)
            throw std::invalid_argument("Duplicate section element");
        cells.push_back(cell);
    }
    std::set<std::pair<std::size_t, std::size_t>> coverage;
    for (std::size_t e = 0; e < source.elements().size(); ++e) {
        std::array<std::size_t, 9> key{};
        for (std::size_t i = 0; i < 8; ++i)
            key[i] = result.source_to_section[source.elements()[e].nodes[bottom[i]]];
        key[8] = block_regions.at(source.element_block_ids()[e]);
        const auto found = sections.find(key);
        if (found == sections.end())
            throw std::invalid_argument("Section connectivity or material changes along z");
        const auto lower =
            std::lower_bound(result.axial.begin(), result.axial.end(), source.nodes()[source.elements()[e].nodes[0]].z);
        const auto layer = static_cast<std::size_t>(lower - result.axial.begin());
        if (layer + 1 >= result.axial.size()
            || result.axial[layer + 1] != source.nodes()[source.elements()[e].nodes[4]].z
            || !coverage.emplace(layer, found->second).second)
            throw std::invalid_argument("Section modal axial mesh overlaps or is nonconforming");
        result.element_section.push_back(found->second);
    }
    if (coverage.size() != (result.axial.size() - 1) * cells.size())
        throw std::invalid_argument("Section modal mesh is missing an extruded cell");
    result.section = std::make_unique<CrossSection>(std::move(nodes), std::move(cells), std::move(regions));
    return result;
}

AxialLayout make_layout(const UnstructuredHex20Mesh& source,
    const Geometry& geometry,
    std::size_t modes,
    const std::vector<ModalEndRegion>& end_regions) {
    const auto nz = geometry.axial.size();
    std::size_t maximum = modes;
    const auto interface = [&](const std::string& name, double fallback) {
        if (name.empty())
            return fallback;
        const auto found = std::find_if(source.node_sets().begin(), source.node_sets().end(), [&](const NodeSet& set) {
            return set.name == name;
        });
        if (found == source.node_sets().end() || found->nodes.empty())
            throw std::invalid_argument("Modal end interface requires a nonempty Exodus node set: " + name);
        const double z = source.nodes()[found->nodes.front()].z;
        std::vector<bool> covered(geometry.section->nodes().size(), false);
        for (auto node : found->nodes) {
            if (source.nodes()[node].z != z)
                throw std::invalid_argument("Modal end interface is not a single axial plane");
            covered[geometry.source_to_section[node]] = true;
        }
        if (std::find(covered.begin(), covered.end(), false) != covered.end()
            || !std::binary_search(geometry.axial.begin(), geometry.axial.end(), z) || z <= geometry.axial.front()
            || z >= geometry.axial.back())
            throw std::invalid_argument("Modal end interface must cover a complete interior axial section");
        return z;
    };
    AxialLayout result;
    result.modes.assign(nz, modes);
    double outer_lower = geometry.axial.front(), outer_upper = geometry.axial.back();
    for (const auto& ends : end_regions) {
        if (ends.mode_count <= modes || ends.mode_count > 768
            || (ends.lower_interface.empty() && ends.upper_interface.empty()))
            throw std::invalid_argument(
                "Local modal count must exceed the full-length count and requires an interface");
        const double lower = interface(ends.lower_interface, geometry.axial.front());
        const double upper = interface(ends.upper_interface, geometry.axial.back());
        outer_lower = std::max(outer_lower, lower);
        outer_upper = std::min(outer_upper, upper);
        if (outer_lower >= outer_upper)
            throw std::invalid_argument("Modal end regions overlap or leave no interior beam interval");
        maximum = std::max(maximum, ends.mode_count);
        for (std::size_t node = 0; node < nz; ++node)
            if (geometry.axial[node] < lower || geometry.axial[node] > upper)
                result.modes[node] = std::max(result.modes[node], ends.mode_count);
    }
    result.indices.resize(3 * maximum * nz, -1);
    for (std::size_t m = 0; m < maximum; ++m)
        for (std::size_t derivative = 0; derivative < 3; ++derivative)
            for (std::size_t node = 0; node < nz; ++node)
                if (m < result.modes[node]) {
                    result.indices[(3 * m + derivative) * nz + node] = static_cast<PetscInt>(result.count++);
                    const bool local = result.modes[node] > modes;
                    result.eliminated.push_back(local);
                    result.retained += local ? 0 : 1;
                }
    return result;
}

std::vector<PetscInt> element_dofs(const Geometry& geometry, std::size_t layer) {
    const auto modes = std::max(geometry.layout.modes[layer], geometry.layout.modes[layer + 1]);
    const auto nz = geometry.axial.size();
    std::vector<PetscInt> result(6 * modes);
    for (std::size_t m = 0; m < modes; ++m)
        for (std::size_t j = 0; j < 6; ++j)
            result[6 * m + j] = geometry.layout.indices[(3 * m + j % 3) * nz + layer + j / 3];
    return result;
}

std::size_t layer_at(const Geometry& geometry, double z) {
    if (z < geometry.axial.front() || z > geometry.axial.back())
        throw std::invalid_argument("Axial sample outside mesh");
    const auto found = std::upper_bound(geometry.axial.begin(), geometry.axial.end(), z);
    return std::min(geometry.axial.size() - 2, static_cast<std::size_t>(found - geometry.axial.begin() - 1));
}

std::vector<double> displacement_row(const Geometry& geometry,
    const ReducedSectionBasis& basis,
    const std::array<std::size_t, 8>& nodes,
    const std::array<double, 8>& shape,
    double z,
    std::size_t component) {
    const auto n = geometry.section->nodes().size();
    const auto layer = layer_at(geometry, z);
    const auto axial = modal_beam_shape(geometry.axial[layer], geometry.axial[layer + 1], z);
    const auto ids = element_dofs(geometry, layer);
    const auto modes = ids.size() / 6;
    std::vector<double> row(geometry.layout.count);
    for (std::size_t m = 0; m < modes; ++m)
        for (std::size_t order = 0; order < 3; ++order) {
            double value = 0.0;
            for (std::size_t i = 0; i < 8; ++i)
                value += shape[i] * basis.modes[m].coefficient[order][component * n + nodes[i]];
            for (std::size_t j = 0; j < 6; ++j)
                if (ids[6 * m + j] >= 0)
                    row[static_cast<std::size_t>(ids[6 * m + j])] += value * axial[order][j];
        }
    return row;
}

std::size_t component(Field field) {
    if (field == Field::displacement_x)
        return 0;
    if (field == Field::displacement_y)
        return 1;
    if (field == Field::displacement_z)
        return 2;
    throw std::invalid_argument("Section modal boundary requires a Cartesian displacement component");
}

struct Constraint final {
    std::vector<double> row;
    double value;
};

void physical_boundaries(const UnstructuredHex20Mesh& source,
    const SpatialDefinition& definition,
    const Geometry& geometry,
    const ReducedSectionBasis& basis,
    std::vector<double>& load,
    std::vector<Constraint>& constraints) {
    std::size_t face_contribution = 0;
    for (const auto& boundary : definition.boundary_conditions) {
        if (!boundary.function.empty() || boundary.use_displaced_geometry)
            throw std::invalid_argument(
                "Section modal boundaries require fixed reference loads without time functions");
        if (boundary.type == BoundaryConditionType::dirichlet) {
            std::set<std::size_t> nodes;
            for (const auto& set : source.node_sets())
                if (set.name == boundary.boundary)
                    nodes.insert(set.nodes.begin(), set.nodes.end());
            if (nodes.empty()) {
                const auto id = source.side_set_block_id(boundary.boundary);
                const auto region = Hex20RegionMesh::from_unstructured_block(source, id);
                for (auto node : region.map_side_set(source, boundary.boundary).displacement_nodes)
                    nodes.insert(region.source_node_ids()[node]);
            }
            if (boundary.field == Field::temperature) {
                for (const auto& region : definition.regions)
                    if (boundary.value != region.initial_temperature)
                        throw std::invalid_argument("Section modal temperature must remain fixed");
                continue;
            }
            for (auto node : nodes) {
                std::array<std::size_t, 8> ids{};
                ids.fill(geometry.source_to_section[node]);
                std::array<double, 8> shape{};
                shape[0] = 1.0;
                constraints.push_back(
                    {displacement_row(geometry, basis, ids, shape, source.nodes()[node].z, component(boundary.field)),
                        boundary.value});
            }
            continue;
        }
        if (boundary.type != BoundaryConditionType::traction && boundary.type != BoundaryConditionType::pressure)
            throw std::invalid_argument(
                "Section modal mechanics supports displacement, traction and pressure boundaries");
        const auto region =
            Hex20RegionMesh::from_unstructured_block(source, source.side_set_block_id(boundary.boundary));
        for (const auto& face : region.map_side_set(source, boundary.boundary).faces) {
            const auto owner = face_contribution++ % static_cast<std::size_t>(PetscGlobalSize);
            if (owner != static_cast<std::size_t>(PetscGlobalRank))
                continue;
            Quad8FaceCoordinates coordinates{};
            std::array<std::size_t, 8> ids{};
            for (std::size_t i = 0; i < 8; ++i) {
                const auto node = region.source_node_ids()[face.nodes[i]];
                coordinates[i] = source.nodes()[node];
                ids[i] = geometry.source_to_section[node];
            }
            for (const auto& qp : make_quad8_face_geometry(coordinates).mechanical_points) {
                double z = 0.0;
                for (std::size_t i = 0; i < 8; ++i)
                    z += qp.displacement_shape[i] * coordinates[i].z;
                if (std::all_of(coordinates.begin(), coordinates.end(), [&](const auto& point) {
                        return point.z == coordinates[0].z;
                    }))
                    z = coordinates[0].z;
                const auto& a = qp.tangent_xi;
                const auto& b = qp.tangent_eta;
                std::array<double, 3> normal{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
                const double area = std::hypot(std::hypot(normal[0], normal[1]), normal[2]);
                for (std::size_t c = 0; c < 3; ++c) {
                    const double force = boundary.type == BoundaryConditionType::pressure
                                             ? -boundary.value * normal[c]
                                             : (c == component(boundary.field) ? boundary.value * area : 0.0);
                    if (force == 0.0)
                        continue;
                    const auto row = displacement_row(geometry, basis, ids, qp.displacement_shape, z, c);
                    for (std::size_t i = 0; i < load.size(); ++i)
                        load[i] += qp.quadrature_weight * force * row[i];
                }
            }
        }
    }
}
} // namespace

SectionModalResult solve_section_modal(const UnstructuredHex20Mesh& source,
    const SpatialDefinition& definition,
    std::size_t full_length_modes,
    const std::vector<ModalEndRegion>& end_regions,
    const SolverOptions& options) {
    const auto start = std::chrono::steady_clock::now();
    const auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    };
    SectionModalResult result;
    PetscBool initialized = PETSC_FALSE;
    check(PetscInitialized(&initialized));
    if (!initialized)
        throw std::logic_error("Section modal solve requires an initialized PETSc session");
    if (options.linear_solver == SolverOptions::LinearSolver::gmres
        || (options.preconditioner != SolverOptions::Preconditioner::automatic
            && options.preconditioner != SolverOptions::Preconditioner::lu)
        || options.mumps_ordering != SolverOptions::MumpsOrdering::automatic || options.maximum_iterations < 1
        || !std::isfinite(options.absolute_tolerance) || options.absolute_tolerance <= 0.0
        || !std::isfinite(options.relative_tolerance) || options.relative_tolerance <= 0.0)
        throw std::invalid_argument("Section modal solve requires direct MUMPS and positive residual tolerances");
    if (full_length_modes < 3 || full_length_modes > 768)
        throw std::invalid_argument("Section modal mechanics requires a mode count in [3,768]");
    // Local regions use a common independent displacement basis. Its classical
    // space includes the two rotations, three Poisson fields and torsion warp;
    // a three-amplitude prefix alone would silently lose those kinematics.
    if (!end_regions.empty() && full_length_modes < 10)
        throw std::invalid_argument("End condensation requires at least 10 full-length modes for the classical space");
    auto geometry = extract(source, definition);
    std::size_t modes = full_length_modes;
    for (const auto& region : end_regions)
        modes = std::max(modes, region.mode_count);
    const auto& section = *geometry.section;
    auto basis = build_reduced_section_basis(section, modes > 3, modes > 4 ? modes - 4 : 0);
    if (modes >= 6)
        basis = make_independent_section_basis(section, make_mixed_bending_basis(std::move(basis)));
    const auto nz = geometry.axial.size();
    if (nz > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()) / (6 * modes))
        throw std::invalid_argument("Section modal axial mesh exceeds PETSc index capacity");
    geometry.layout = make_layout(source, geometry, full_length_modes, end_regions);
    std::map<std::size_t, ReducedSectionBasis> element_bases;
    std::map<std::size_t, ModalSectionKinematics> section_kinematics;
    std::map<std::size_t, ModalSectionLinearization> section_linearizations;
    for (auto local_count : geometry.layout.modes)
        if (element_bases.count(local_count) == 0) {
            auto prefix = basis;
            prefix.modes.resize(local_count);
            section_kinematics.emplace(local_count, sample_modal_section(section, prefix));
            section_linearizations.emplace(local_count, linearize_modal_section(section, prefix));
            element_bases.emplace(local_count, std::move(prefix));
        }
    const auto count = geometry.layout.count;
    const auto n = static_cast<PetscInt>(count);
    const int rank = PetscGlobalRank, ranks = PetscGlobalSize;
    Objects objects;
    PetscInt local_rows = PETSC_DECIDE, global_rows = n;
    check(PetscSplitOwnership(PETSC_COMM_WORLD, &local_rows, &global_rows));
    // Each axial node couples to all three jets of every mode at itself and
    // at most two neighbors. Preallocate this known graph before inserting
    // element matrices; zero preallocation otherwise reallocates the entire
    // growing AIJ storage repeatedly on a long production extrusion.
    const auto row_entries = static_cast<PetscInt>(9 * modes);
    check(MatCreateAIJ(PETSC_COMM_WORLD,
        local_rows,
        local_rows,
        n,
        n,
        std::min(local_rows, row_entries),
        nullptr,
        std::min(n - local_rows, row_entries),
        nullptr,
        &objects.stiffness));
    check(MatSetOption(objects.stiffness, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    result.section_preprocessing_seconds = elapsed();
    for (std::size_t e = 0; e + 1 < nz; ++e) {
        if (e % static_cast<std::size_t>(ranks) != static_cast<std::size_t>(rank))
            continue;
        const auto element = make_modal_beam_element(geometry.axial[e], geometry.axial[e + 1]);
        const auto ids = element_dofs(geometry, e);
        const auto local = modal_beam_linear_stiffness(section_linearizations.at(ids.size() / 6), element);
        const auto size = static_cast<PetscInt>(ids.size());
        check(MatSetValues(objects.stiffness, size, ids.data(), size, ids.data(), local.data(), ADD_VALUES));
    }
    check(MatAssemblyBegin(objects.stiffness, MAT_FINAL_ASSEMBLY));
    check(MatAssemblyEnd(objects.stiffness, MAT_FINAL_ASSEMBLY));
    check(MatCreateVecs(objects.stiffness, &objects.diagonal, nullptr));
    check(MatGetDiagonal(objects.stiffness, objects.diagonal));
    auto scaling = gather(objects.diagonal);
    for (double& value : scaling) {
        if (!(value > 0.0) || !std::isfinite(value))
            throw std::runtime_error("Modal stiffness has an inactive or invalid unknown");
        value = 1.0 / std::sqrt(value);
    }
    std::vector<double> load(count);
    std::vector<Constraint> physical, orthogonal;
    physical_boundaries(source, definition, geometry, basis, load, physical);
    // Each physical face, like each axial element, is evaluated on one rank.
    // Only the small generalized load vector is gathered for the constrained solve.
    check(MatCreateVecs(objects.stiffness, &objects.physical_load, nullptr));
    check(VecSet(objects.physical_load, 0.0));
    for (std::size_t i = 0; i < count; ++i)
        if (load[i] != 0.0)
            check(VecSetValue(objects.physical_load, static_cast<PetscInt>(i), load[i], ADD_VALUES));
    check(VecAssemblyBegin(objects.physical_load));
    check(VecAssemblyEnd(objects.physical_load));
    load = gather(objects.physical_load);
    // Exact physical displacement constraints, not an artificial clamp on modal derivatives.
    // Diagonal equilibration followed by row-pivoted, twice-reorthogonalized QR.
    // Pivoting matters: accepting a nearly dependent early row can amplify tiny
    // section auxiliary-solve errors and invent additional end constraints.
    std::vector<Constraint> candidates;
    std::vector<double> physical_scale;
    std::vector<std::size_t> active_columns;
    for (std::size_t i = 0; i < count; ++i)
        if (std::any_of(physical.begin(), physical.end(), [i](const Constraint& constraint) {
                return constraint.row[i] != 0.0;
            }))
            active_columns.push_back(i);
    // End constraints touch only end amplitudes. Keep every nonzero column,
    // however small, but avoid sweeping thousands of identically zero interior
    // columns during each QR pass. Expand the result back to the global layout
    // before assembly. This is an exact graph compression, not a rank tolerance.
    for (const auto& constraint : physical) {
        std::vector<double> row(active_columns.size());
        for (std::size_t i = 0; i < row.size(); ++i)
            row[i] = constraint.row[active_columns[i]] * scaling[active_columns[i]];
        const double norm = std::sqrt(std::inner_product(row.begin(), row.end(), row.begin(), 0.0));
        if (norm == 0.0) {
            if (constraint.value != 0.0)
                throw std::invalid_argument("Physical displacement cannot be represented by the selected modes");
            continue;
        }
        for (double& value : row)
            value /= norm;
        candidates.push_back({std::move(row), constraint.value / norm});
        physical_scale.push_back(norm);
    }
    while (!candidates.empty()) {
        double largest = 0.0;
        std::size_t pivot = 0;
        Constraint selected;
        for (std::size_t candidate = 0; candidate < candidates.size(); ++candidate) {
            auto row = candidates[candidate].row;
            double value = candidates[candidate].value;
            for (int pass = 0; pass < 2; ++pass)
                for (const auto& previous : orthogonal) {
                    const double projection = std::inner_product(row.begin(), row.end(), previous.row.begin(), 0.0);
                    for (std::size_t i = 0; i < row.size(); ++i)
                        row[i] -= projection * previous.row[i];
                    value -= projection * previous.value;
                }
            const double norm = std::sqrt(std::inner_product(row.begin(), row.end(), row.begin(), 0.0));
            if (norm < 1.0e-10 && std::abs(value) * physical_scale[candidate] > 1.0e-10)
                throw std::invalid_argument("Inconsistent projected displacement constraints");
            if (norm > largest) {
                largest = norm;
                pivot = candidate;
                selected = {std::move(row), value};
            }
        }
        if (largest < 1.0e-10)
            break;
        for (double& entry : selected.row)
            entry /= largest;
        selected.value /= largest;
        orthogonal.push_back(std::move(selected));
        candidates.erase(candidates.begin() + static_cast<std::ptrdiff_t>(pivot));
        physical_scale.erase(physical_scale.begin() + static_cast<std::ptrdiff_t>(pivot));
    }
    for (auto& constraint : orthogonal) {
        std::vector<double> row(count);
        for (std::size_t i = 0; i < active_columns.size(); ++i)
            row[active_columns[i]] = constraint.row[i];
        constraint.row = std::move(row);
    }
    const auto total = static_cast<PetscInt>(count + orthogonal.size());
    PetscInt local_total = PETSC_DECIDE, global_total = total;
    check(PetscSplitOwnership(PETSC_COMM_WORLD, &local_total, &global_total));
    check(MatCreateAIJ(PETSC_COMM_WORLD,
        local_total,
        local_total,
        total,
        total,
        std::min(local_total, row_entries),
        nullptr,
        std::min(total - local_total, row_entries),
        nullptr,
        &objects.constrained));
    check(MatSetOption(objects.constrained, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    PetscInt begin, end;
    check(MatGetOwnershipRange(objects.stiffness, &begin, &end));
    for (PetscInt i = begin; i < end; ++i) {
        PetscInt size;
        const PetscInt* columns;
        const PetscScalar* entries;
        check(MatGetRow(objects.stiffness, i, &size, &columns, &entries));
        std::vector<double> values(static_cast<std::size_t>(size));
        for (PetscInt j = 0; j < size; ++j)
            values[static_cast<std::size_t>(j)] =
                entries[j] * scaling[static_cast<std::size_t>(i)] * scaling[static_cast<std::size_t>(columns[j])];
        check(MatSetValues(objects.constrained, 1, &i, size, columns, values.data(), INSERT_VALUES));
        check(MatRestoreRow(objects.stiffness, i, &size, &columns, &entries));
    }
    if (rank == 0)
        for (std::size_t c = 0; c < orthogonal.size(); ++c) {
            const auto row = static_cast<PetscInt>(count + c);
            for (std::size_t i = 0; i < count; ++i)
                if (orthogonal[c].row[i] != 0.0) {
                    const auto column = static_cast<PetscInt>(i);
                    const double value = orthogonal[c].row[i];
                    check(MatSetValue(objects.constrained, row, column, value, INSERT_VALUES));
                    check(MatSetValue(objects.constrained, column, row, value, INSERT_VALUES));
                }
            check(MatSetValue(objects.constrained, row, row, 0.0, INSERT_VALUES));
        }
    check(MatAssemblyBegin(objects.constrained, MAT_FINAL_ASSEMBLY));
    check(MatAssemblyEnd(objects.constrained, MAT_FINAL_ASSEMBLY));
    // The scaled saddle matrix now owns every stiffness entry. Release the
    // unscaled copy before factorization; retaining both doubles the dominant
    // assembly storage on sections with many independent amplitudes.
    check(MatDestroy(&objects.stiffness));
    check(VecDestroy(&objects.diagonal));
    check(VecDestroy(&objects.physical_load));
    check(MatCreateVecs(objects.constrained, &objects.solution, &objects.rhs));
    if (rank == 0) {
        for (std::size_t i = 0; i < count; ++i)
            check(VecSetValue(objects.rhs, static_cast<PetscInt>(i), load[i] * scaling[i], INSERT_VALUES));
        for (std::size_t c = 0; c < orthogonal.size(); ++c)
            check(VecSetValue(objects.rhs, static_cast<PetscInt>(count + c), orthogonal[c].value, INSERT_VALUES));
    }
    check(VecAssemblyBegin(objects.rhs));
    check(VecAssemblyEnd(objects.rhs));
    std::unique_ptr<ModalStaticCondensation> condensed;
    result.global_constraint_count = orthogonal.size();
    result.assembly_seconds = elapsed() - result.section_preprocessing_seconds;
    const double condensation_start = elapsed();
    if (!end_regions.empty()) {
        auto eliminated = geometry.layout.eliminated;
        for (const auto& constraint : orthogonal) {
            bool local = false, retained = false;
            for (std::size_t i = 0; i < count; ++i)
                if (constraint.row[i] != 0.0) {
                    local = local || eliminated[i];
                    retained = retained || !eliminated[i];
                }
            if (local && retained)
                throw std::invalid_argument("Physical displacement constraint crosses a modal condensation interface");
            eliminated.push_back(local);
            if (local)
                --result.global_constraint_count;
        }
        condensed = std::make_unique<ModalStaticCondensation>(objects.constrained, eliminated);
    } else {
        check(KSPCreate(PETSC_COMM_WORLD, &objects.solver));
        check(KSPSetOperators(objects.solver, objects.constrained, objects.constrained));
        check(KSPSetType(objects.solver, KSPPREONLY));
        PC pc;
        check(KSPGetPC(objects.solver, &pc));
        check(PCSetType(pc, PCLU));
        check(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    }
    result.condensation_seconds = elapsed() - condensation_start;
    const double solve_start = elapsed();
    const auto solve = [&](Vec rhs, Vec solution) {
        if (condensed)
            condensed->solve(rhs, solution);
        else {
            check(KSPSolve(objects.solver, rhs, solution));
            KSPConvergedReason reason;
            check(KSPGetConvergedReason(objects.solver, &reason));
            if (reason <= 0)
                throw std::runtime_error("Modal axial solve failed; check physical rigid-body constraints");
        }
    };
    solve(objects.rhs, objects.solution);
    check(VecDuplicate(objects.rhs, &objects.defect));
    // Inactive interface jets use negative indices, as for homogeneous
    // constraints. PETSc vectors require this explicit option, unlike matrices.
    check(VecSetOption(objects.defect, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE));
    check(VecDuplicate(objects.rhs, &objects.correction));
    PetscReal defect, norm;
    check(VecNorm(objects.rhs, NORM_2, &norm));
    // On a refined slender beam K*x subtracts very large entries to recover
    // small bending forces. Refine the direct solution using the original
    // material-point virtual work, not a repeatedly rounded matrix product.
    // The same LU factors solve each correction. Element residuals, like their
    // tangents, are evaluated by exactly one MPI rank. This is linear defect
    // correction, not a change of stiffness, load, or acceptance tolerance.
    for (std::size_t iteration = 0; iteration <= 6; ++iteration) {
        const auto current = gather(objects.solution);
        check(VecSet(objects.defect, 0.0));
        for (std::size_t e = 0; e + 1 < nz; ++e) {
            if (e % static_cast<std::size_t>(ranks) != static_cast<std::size_t>(rank))
                continue;
            const auto ids = element_dofs(geometry, e);
            std::vector<double> local_state(ids.size());
            for (std::size_t i = 0; i < ids.size(); ++i) {
                if (ids[i] < 0)
                    continue;
                const auto global = static_cast<std::size_t>(ids[i]);
                local_state[i] = current[global] * scaling[global];
            }
            auto response = evaluate_modal_beam(section,
                element_bases.at(ids.size() / 6),
                make_modal_beam_element(geometry.axial[e], geometry.axial[e + 1]),
                local_state,
                false,
                &section_kinematics.at(ids.size() / 6));
            for (std::size_t i = 0; i < ids.size(); ++i)
                if (ids[i] >= 0)
                    response.residual[i] *= scaling[static_cast<std::size_t>(ids[i])];
            check(VecSetValues(objects.defect,
                static_cast<PetscInt>(ids.size()),
                ids.data(),
                response.residual.data(),
                ADD_VALUES));
        }
        if (rank == 0) {
            for (std::size_t i = 0; i < count; ++i) {
                double value = -load[i] * scaling[i];
                for (std::size_t c = 0; c < orthogonal.size(); ++c)
                    value += orthogonal[c].row[i] * current[count + c];
                check(VecSetValue(objects.defect, static_cast<PetscInt>(i), value, ADD_VALUES));
            }
            for (std::size_t c = 0; c < orthogonal.size(); ++c) {
                const double value =
                    std::inner_product(orthogonal[c].row.begin(), orthogonal[c].row.end(), current.begin(), 0.0)
                    - orthogonal[c].value;
                check(VecSetValue(objects.defect, static_cast<PetscInt>(count + c), value, ADD_VALUES));
            }
        }
        check(VecAssemblyBegin(objects.defect));
        check(VecAssemblyEnd(objects.defect));
        check(VecNorm(objects.defect, NORM_2, &defect));
        PetscScalar work_defect;
        check(VecDot(objects.solution, objects.defect, &work_defect));
        double external_work = 0.0;
        for (std::size_t i = 0; i < count; ++i)
            external_work += load[i] * scaling[i] * current[i];
        for (std::size_t c = 0; c < orthogonal.size(); ++c)
            external_work -= orthogonal[c].value * current[count + c];
        result.refinement_iterations = iteration;
        if ((defect <= std::min(options.absolute_tolerance, 1.0e-9 * norm)
                && std::abs(work_defect) <= 1.0e-10 * std::abs(external_work))
            || iteration == 6)
            break;
        solve(objects.defect, objects.correction);
        check(VecAXPY(objects.solution, -1.0, objects.correction));
    }
    check(MatMult(objects.constrained, objects.solution, objects.correction));
    check(VecAXPY(objects.correction, -1.0, objects.rhs));
    PetscReal algebraic_defect;
    check(VecNorm(objects.correction, NORM_2, &algebraic_defect));
    result.algebraic_relative_residual = norm > 0.0 ? algebraic_defect / norm : algebraic_defect;
    result.mode_count = full_length_modes;
    result.basis_mode_count = modes;
    result.global_dof_count = geometry.layout.retained;
    result.condensed_dof_count = count - geometry.layout.retained;
    result.nodal_mode_counts = geometry.layout.modes;
    result.amplitude_map.resize(geometry.layout.indices.size(), count);
    for (std::size_t i = 0; i < geometry.layout.indices.size(); ++i)
        if (geometry.layout.indices[i] >= 0)
            result.amplitude_map[i] = static_cast<std::size_t>(geometry.layout.indices[i]);
    for (const auto& mode : basis.modes) {
        result.poisson_modes += mode.kind == SectionMode::Kind::poisson_relaxation ? 1 : 0;
        result.distortion_modes += mode.kind == SectionMode::Kind::distortion ? 1 : 0;
        result.axial_warping_modes += mode.kind == SectionMode::Kind::axial_warping ? 1 : 0;
        result.shear_free_modes += mode.kind == SectionMode::Kind::shear_free_distortion ? 1 : 0;
        result.transverse_corrector_modes += mode.kind == SectionMode::Kind::transverse_corrector ? 1 : 0;
        result.shear_modes +=
            (mode.kind == SectionMode::Kind::shear_x || mode.kind == SectionMode::Kind::shear_y) ? 1 : 0;
    }
    result.axial_nodes = nz;
    result.axial_coordinates = geometry.axial;
    result.constraint_count = orthogonal.size();
    result.equilibrium_relative_residual = norm > 0.0 ? defect / norm : defect;
    if (!std::isfinite(defect) || !std::isfinite(norm)
        || defect > std::max(options.absolute_tolerance, options.relative_tolerance * norm))
        throw std::runtime_error("Modal axial equilibrium residual exceeds tolerance");
    const auto solution = gather(objects.solution);
    result.amplitudes.resize(count);
    for (std::size_t i = 0; i < count; ++i)
        result.amplitudes[i] = solution[i] * scaling[i];
    result.external_work = std::inner_product(load.begin(), load.end(), result.amplitudes.begin(), 0.0);
    for (std::size_t c = 0; c < orthogonal.size(); ++c)
        result.external_work -= orthogonal[c].value * solution[count + c];
    for (const auto& constraint : physical) {
        const double actual =
            std::inner_product(constraint.row.begin(), constraint.row.end(), result.amplitudes.begin(), 0.0);
        result.constraint_maximum_error =
            std::max(result.constraint_maximum_error, std::abs(actual - constraint.value));
    }
    if (result.constraint_maximum_error > 1.0e-9)
        throw std::runtime_error("Reconstructed physical displacement constraint failed");
    result.solve_and_refinement_seconds = elapsed() - solve_start;
    const double recovery_start = elapsed();
    // Replicated postprocessing is read-only; each stiffness contribution above has exactly one rank owner.
    for (std::size_t e = 0; e + 1 < nz; ++e) {
        const auto element = make_modal_beam_element(geometry.axial[e], geometry.axial[e + 1]);
        const auto ids = element_dofs(geometry, e);
        std::vector<double> state(ids.size());
        for (std::size_t i = 0; i < ids.size(); ++i)
            if (ids[i] >= 0)
                state[i] = result.amplitudes[static_cast<std::size_t>(ids[i])];
        const auto response = evaluate_modal_beam(section,
            element_bases.at(ids.size() / 6),
            element,
            state,
            false,
            &section_kinematics.at(ids.size() / 6));
        result.energy += response.energy;
        for (std::size_t z = 0; z < element.points.size(); ++z) {
            SectionModalResultant resultant;
            resultant.z = element.points[z].z;
            for (std::size_t q = 0; q < section.points().size(); ++q) {
                const auto& point = section.points()[q];
                auto position = point.position;
                position.z = element.points[z].z;
                const auto index = z * section.points().size() + q;
                result.samples.push_back({e,
                    q,
                    position,
                    point.weight * element.points[z].weight,
                    response.strain[index],
                    response.stress[index],
                    response.displacement[index]});
                const auto& stress = response.stress[index];
                const double x = point.position.x - section.elastic_center().x;
                const double y = point.position.y - section.elastic_center().y;
                resultant.axial_force += point.weight * stress[2];
                resultant.moment_x += point.weight * y * stress[2];
                resultant.moment_y -= point.weight * x * stress[2];
                resultant.torque += point.weight * (x * stress[4] - y * stress[5]);
                for (std::size_t c = 0; c < 6; ++c)
                    resultant.energy_per_length +=
                        0.5 * point.weight * (c < 3 ? 1.0 : 2.0) * stress[c] * response.strain[index][c];
            }
            result.resultants.push_back(resultant);
        }
    }
    for (std::size_t node = 0; node < source.nodes().size(); ++node) {
        std::array<std::size_t, 8> ids{};
        ids.fill(geometry.source_to_section[node]);
        std::array<double, 8> shape{};
        shape[0] = 1.0;
        std::array<double, 3> u{};
        for (std::size_t c = 0; c < 3; ++c) {
            const auto row = displacement_row(geometry, basis, ids, shape, source.nodes()[node].z, c);
            u[c] = std::inner_product(row.begin(), row.end(), result.amplitudes.begin(), 0.0);
        }
        result.displacement.push_back({u[0], u[1], u[2]});
    }
    result.field_recovery_seconds = elapsed() - recovery_start;
    return result;
}
} // namespace fuelsim
