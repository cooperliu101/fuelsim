#include "cartesian3d_assembly.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

namespace fuelsim::cartesian {
namespace {
CartesianPoint3 face_centroid(const Quad4FaceCoordinates& coordinates) {
    CartesianPoint3 result{0.0, 0.0, 0.0};
    for (const CartesianPoint3& point : coordinates) {
        result.x += 0.25 * point.x;
        result.y += 0.25 * point.y;
        result.z += 0.25 * point.z;
    }
    return result;
}

CartesianPoint3 element_centroid(const Hex8RegionMesh& mesh, std::size_t element_index) {
    CartesianPoint3 result{0.0, 0.0, 0.0};
    for (std::size_t node : mesh.elements().at(element_index).nodes) {
        result.x += 0.125 * mesh.nodes().at(node).x;
        result.y += 0.125 * mesh.nodes().at(node).y;
        result.z += 0.125 * mesh.nodes().at(node).z;
    }
    return result;
}

Quad4FaceCoordinates face_coordinates(const Hex8RegionMesh& mesh, const Quad4FaceElement& face) {
    Quad4FaceCoordinates result{};
    for (std::size_t node = 0; node < 4; ++node) result[node] = mesh.nodes().at(face.nodes[node]);
    return result;
}

CartesianPoint3 subtract(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

CartesianPoint3 cross(const CartesianPoint3& first, const CartesianPoint3& second) {
    return {first.y * second.z - first.z * second.y, first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

double dot(const CartesianPoint3& first, const CartesianPoint3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

using Matrix8 = std::array<double, 64>;
using Vector8 = std::array<double, 8>;

const std::array<std::array<double, 2>, 8>& abaqus_quad8_constraint_locations() {
    // H20.28 identifies the effective centers of Abaqus/Standard's default
    // quadratic secondary constraint regions.  Coordinates are in the Q8
    // parent domain and follow the face-node order.
    static const std::array<std::array<double, 2>, 8> value = {
        {{-0.75, -0.75}, {0.75, -0.75}, {0.75, 0.75}, {-0.75, 0.75}, {0.0, -0.5}, {0.5, 0.0}, {0.0, 0.5}, {-0.5, 0.0}}};
    return value;
}

const Matrix8& abaqus_quad8_averaging() {
    // H20.26 and H20.28 identify Abaqus/Standard's default quadratic
    // small-sliding C3D20 secondary-face averaging operator.  This is an
    // observed compatibility operator, not a classical mortar dual basis.
    static const Matrix8 value = {{
        0.4051580465416001,
        -0.08720743312758768,
        -0.02922569057450125,
        -0.08720743312758091,
        0.3574287590139803,
        0.0418124961300837,
        0.0418124961300837,
        0.3574287590139803,
        -0.08720743312758768,
        0.4051580465416001,
        -0.08720743312758768,
        -0.02922569057450803,
        0.3574287590139803,
        0.3574287590139803,
        0.0418124961300837,
        0.04181249613007693,
        -0.0292256905745148,
        -0.08720743312759446,
        0.4051580465415933,
        -0.08720743312759446,
        0.04181249613007693,
        0.3574287590139803,
        0.3574287590139803,
        0.04181249613007693,
        -0.08720743312758091,
        -0.02922569057450125,
        -0.08720743312758091,
        0.4051580465415933,
        0.04181249613007693,
        0.04181249613007693,
        0.3574287590139803,
        0.3574287590139803,
        -0.08706059524915631,
        -0.08706059524915631,
        -0.1330911537220558,
        -0.1330911537220558,
        0.6107531349427593,
        0.313368055555558,
        0.2028142518885394,
        0.313368055555558,
        -0.1330911537220219,
        -0.08706059524912242,
        -0.08706059524912242,
        -0.1330911537220287,
        0.3133680555555783,
        0.6107531349427796,
        0.3133680555555716,
        0.2028142518885123,
        -0.1330911537220355,
        -0.1330911537220355,
        -0.08706059524913598,
        -0.08706059524913598,
        0.2028142518885258,
        0.3133680555555377,
        0.610753134942739,
        0.3133680555555445,
        -0.08706059524916308,
        -0.1330911537220558,
        -0.1330911537220558,
        -0.08706059524916308,
        0.3133680555555648,
        0.2028142518885462,
        0.313368055555558,
        0.6107531349427593,
    }};
    return value;
}

Vector8 solve8(Matrix8 matrix, Vector8 right_hand_side) {
    for (std::size_t column = 0; column < 8; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < 8; ++row)
            if (std::abs(matrix[row * 8 + column]) > std::abs(matrix[pivot * 8 + column])) pivot = row;
        if (!(std::abs(matrix[pivot * 8 + column]) > 1.0e-14))
            throw std::invalid_argument("HEX20 averaged contact has a singular secondary face Gram matrix");
        if (pivot != column) {
            for (std::size_t entry = column; entry < 8; ++entry)
                std::swap(matrix[column * 8 + entry], matrix[pivot * 8 + entry]);
            std::swap(right_hand_side[column], right_hand_side[pivot]);
        }
        const double diagonal = matrix[column * 8 + column];
        for (std::size_t row = column + 1; row < 8; ++row) {
            const double factor = matrix[row * 8 + column] / diagonal;
            for (std::size_t entry = column; entry < 8; ++entry)
                matrix[row * 8 + entry] -= factor * matrix[column * 8 + entry];
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    Vector8 result{};
    for (std::size_t reverse = 0; reverse < 8; ++reverse) {
        const std::size_t row = 7 - reverse;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < 8; ++column) value -= matrix[row * 8 + column] * result[column];
        result[row] = value / matrix[row * 8 + row];
    }
    return result;
}

double reference_measure(const Quad8FaceMechanicalQuadraturePoint& point) {
    const CartesianPoint3 area_vector = cross(point.tangent_xi, point.tangent_eta);
    return std::sqrt(dot(area_vector, area_vector));
}

std::array<double, 3> thermal_search_point(
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    std::array<double, 3> result{};
    for (std::size_t node = 0; node < 4; ++node) {
        result[0] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].x + state[8 + node]);
        result[1] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].y + state[16 + node]);
        result[2] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].z + state[24 + node]);
    }
    return result;
}

std::array<double, 3> mechanical_search_point(
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    const std::size_t node = geometry.secondary_local_node;
    return {geometry.secondary_coordinates[node].x + state[8 + node],
        geometry.secondary_coordinates[node].y + state[16 + node],
        geometry.secondary_coordinates[node].z + state[24 + node]};
}

std::array<double, 3> hex20_thermal_search_point(
    const Quad8ToQuad8HeatGeometry& geometry, const Quad8SurfaceContactLocalValues& state) {
    std::array<double, 3> result{};
    for (std::size_t node = 0; node < 8; ++node) {
        result[0] +=
            geometry.secondary_displacement_shape[node] * (geometry.secondary_coordinates[node].x + state[8 + node]);
        result[1] +=
            geometry.secondary_displacement_shape[node] * (geometry.secondary_coordinates[node].y + state[24 + node]);
        result[2] +=
            geometry.secondary_displacement_shape[node] * (geometry.secondary_coordinates[node].z + state[40 + node]);
    }
    return result;
}

std::array<double, 3> hex20_mechanical_search_point(
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalValues& state) {
    const std::size_t node = geometry.secondary_local_node;
    return {geometry.secondary_coordinates[node].x + state[8 + node],
        geometry.secondary_coordinates[node].y + state[24 + node],
        geometry.secondary_coordinates[node].z + state[40 + node]};
}

std::array<double, 3> hex20_mechanical_search_point(
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalValues& state) {
    std::array<double, 3> result{};
    for (std::size_t node = 0; node < 8; ++node) {
        result[0] +=
            geometry.secondary_displacement_shape[node] * (geometry.secondary_coordinates[node].x + state[8 + node]);
        result[1] +=
            geometry.secondary_displacement_shape[node] * (geometry.secondary_coordinates[node].y + state[24 + node]);
        result[2] +=
            geometry.secondary_displacement_shape[node] * (geometry.secondary_coordinates[node].z + state[40 + node]);
    }
    return result;
}

double normal_orientation(const Quad4FaceCoordinates& primary_coordinates,
    const CartesianPoint3& secondary_parent_centroid, const CartesianPoint3& primary_parent_centroid) {
    const CartesianPoint3 tangent_xi = subtract(primary_coordinates[1], primary_coordinates[0]),
                          tangent_eta = subtract(primary_coordinates[3], primary_coordinates[0]),
                          area = cross(tangent_xi, tangent_eta),
                          material_direction = subtract(primary_parent_centroid, secondary_parent_centroid);
    const double measure = std::sqrt(dot(area, area)), orientation = dot(area, material_direction);
    if (!std::isfinite(measure) || !(measure > 0.0) || !std::isfinite(orientation) || orientation == 0.0)
        throw std::invalid_argument(
            "Three-dimensional contact faces require nondegenerate opposing material-side centroids");
    return orientation > 0.0 ? 1.0 : -1.0;
}

double minimum_normal_length(const Hex8RegionMesh& mesh, const Hex8RegionBoundary& boundary, const std::string& name) {
    double result = std::numeric_limits<double>::infinity();
    for (const Quad4FaceElement& face : boundary.faces) {
        const Quad4FaceCoordinates coordinates = face_coordinates(mesh, face);
        const CartesianPoint3 face_center = face_centroid(coordinates),
                              parent = element_centroid(mesh, face.parent_element),
                              first = subtract(coordinates[1], coordinates[0]),
                              second = subtract(coordinates[3], coordinates[0]), area = cross(first, second);
        const double measure = std::sqrt(dot(area, area));
        if (!std::isfinite(measure) || !(measure > 0.0))
            throw std::invalid_argument("Three-dimensional contact boundary has a degenerate face: " + name);
        const double height = 2.0 * std::abs(dot(subtract(parent, face_center), area)) / measure;
        if (!std::isfinite(height) || !(height > 0.0))
            throw std::invalid_argument("Three-dimensional contact boundary has a nonpositive normal length: " + name);
        result = std::min(result, height);
    }
    if (!std::isfinite(result)) throw std::invalid_argument("Three-dimensional contact boundary is empty: " + name);
    return result;
}

CartesianPoint3 hex20_element_centroid(const Hex20RegionMesh& mesh, std::size_t element_index) {
    CartesianPoint3 result{0.0, 0.0, 0.0};
    const Hex20Element& element = mesh.elements().at(element_index);
    for (std::size_t node = 0; node < 8; ++node) {
        result.x += 0.125 * mesh.nodes().at(element.nodes[node]).x;
        result.y += 0.125 * mesh.nodes().at(element.nodes[node]).y;
        result.z += 0.125 * mesh.nodes().at(element.nodes[node]).z;
    }
    return result;
}

Quad8FaceCoordinates face_coordinates(const Hex20RegionMesh& mesh, const Quad8FaceElement& face) {
    Quad8FaceCoordinates result{};
    for (std::size_t node = 0; node < 8; ++node) result[node] = mesh.nodes().at(face.nodes[node]);
    return result;
}

CartesianPoint3 quad8_face_centroid(const Quad8FaceCoordinates& coordinates) {
    CartesianPoint3 result{0.0, 0.0, 0.0};
    for (std::size_t node = 0; node < 4; ++node) {
        result.x -= 0.25 * coordinates[node].x;
        result.y -= 0.25 * coordinates[node].y;
        result.z -= 0.25 * coordinates[node].z;
    }
    for (std::size_t node = 4; node < 8; ++node) {
        result.x += 0.5 * coordinates[node].x;
        result.y += 0.5 * coordinates[node].y;
        result.z += 0.5 * coordinates[node].z;
    }
    return result;
}

double normal_orientation(const Quad8FaceCoordinates& primary_coordinates,
    const CartesianPoint3& secondary_parent_centroid, const CartesianPoint3& primary_parent_centroid) {
    const Quad8FaceGeometry geometry = make_quad8_face_geometry(primary_coordinates);
    const Quad8FaceMechanicalQuadraturePoint& center = geometry.mechanical_points[4];
    const CartesianPoint3 area{center.tangent_xi.y * center.tangent_eta.z - center.tangent_xi.z * center.tangent_eta.y,
        center.tangent_xi.z * center.tangent_eta.x - center.tangent_xi.x * center.tangent_eta.z,
        center.tangent_xi.x * center.tangent_eta.y - center.tangent_xi.y * center.tangent_eta.x};
    const CartesianPoint3 material_direction = subtract(primary_parent_centroid, secondary_parent_centroid);
    const double measure = std::sqrt(dot(area, area)), orientation = dot(area, material_direction);
    if (!std::isfinite(measure) || !(measure > 0.0) || !std::isfinite(orientation) || orientation == 0.0)
        throw std::invalid_argument("HEX20 contact faces require nondegenerate opposing material-side centroids");
    return orientation > 0.0 ? 1.0 : -1.0;
}

double minimum_normal_length(
    const Hex20RegionMesh& mesh, const Hex20RegionBoundary& boundary, const std::string& name) {
    double result = std::numeric_limits<double>::infinity();
    for (const Quad8FaceElement& face : boundary.faces) {
        const Quad8FaceCoordinates coordinates = face_coordinates(mesh, face);
        const CartesianPoint3 face_center = quad8_face_centroid(coordinates),
                              parent = hex20_element_centroid(mesh, face.parent_element);
        const Quad8FaceGeometry geometry = make_quad8_face_geometry(coordinates);
        const Quad8FaceMechanicalQuadraturePoint& center = geometry.mechanical_points[4];
        const CartesianPoint3 area{
            center.tangent_xi.y * center.tangent_eta.z - center.tangent_xi.z * center.tangent_eta.y,
            center.tangent_xi.z * center.tangent_eta.x - center.tangent_xi.x * center.tangent_eta.z,
            center.tangent_xi.x * center.tangent_eta.y - center.tangent_xi.y * center.tangent_eta.x};
        const double measure = std::sqrt(dot(area, area));
        if (!std::isfinite(measure) || !(measure > 0.0))
            throw std::invalid_argument("HEX20 contact boundary has a degenerate face: " + name);
        const double height = 2.0 * std::abs(dot(subtract(parent, face_center), area)) / measure;
        if (!std::isfinite(height) || !(height > 0.0))
            throw std::invalid_argument("HEX20 contact boundary has a nonpositive normal length: " + name);
        result = std::min(result, height);
    }
    if (!std::isfinite(result)) throw std::invalid_argument("HEX20 contact boundary is empty: " + name);
    return result;
}

spatial_detail::ContactSearchBox quad8_search_box(const Quad8FaceCoordinates& coordinates,
    const std::vector<double>& state, const std::array<std::array<std::size_t, 3>, 8>& displacement_dofs) {
    std::array<CartesianPoint3, 8> current = coordinates;
    for (std::size_t node = 0; node < 8; ++node) {
        current[node].x += state[displacement_dofs[node][0]];
        current[node].y += state[displacement_dofs[node][1]];
        current[node].z += state[displacement_dofs[node][2]];
    }
    spatial_detail::ContactSearchBox box;
    box.minimum.fill(std::numeric_limits<double>::infinity());
    box.maximum.fill(-std::numeric_limits<double>::infinity());
    for (std::size_t component = 0; component < 3; ++component) {
        const auto coordinate = [component](const CartesianPoint3& point) {
            return component == 0 ? point.x : (component == 1 ? point.y : point.z);
        };
        std::array<std::array<double, 3>, 3> lagrange{};
        lagrange[0][0] = coordinate(current[0]);
        lagrange[2][0] = coordinate(current[1]);
        lagrange[2][2] = coordinate(current[2]);
        lagrange[0][2] = coordinate(current[3]);
        lagrange[1][0] = coordinate(current[4]);
        lagrange[2][1] = coordinate(current[5]);
        lagrange[1][2] = coordinate(current[6]);
        lagrange[0][1] = coordinate(current[7]);
        lagrange[1][1] = coordinate(quad8_face_centroid(current));
        std::array<std::array<double, 3>, 3> x_bernstein{};
        for (std::size_t j = 0; j < 3; ++j) {
            x_bernstein[0][j] = lagrange[0][j];
            x_bernstein[2][j] = lagrange[2][j];
            x_bernstein[1][j] = 2.0 * lagrange[1][j] - 0.5 * (lagrange[0][j] + lagrange[2][j]);
        }
        for (std::size_t i = 0; i < 3; ++i) {
            const double values[3] = {x_bernstein[i][0],
                2.0 * x_bernstein[i][1] - 0.5 * (x_bernstein[i][0] + x_bernstein[i][2]), x_bernstein[i][2]};
            for (double value : values) {
                box.minimum[component] = std::min(box.minimum[component], value);
                box.maximum[component] = std::max(box.maximum[component], value);
            }
        }
    }
    return box;
}

void set_pattern_block(std::vector<unsigned char>& pattern, std::size_t size, std::size_t row_begin,
    std::size_t row_end, std::size_t column_begin, std::size_t column_end) {
    for (std::size_t row = row_begin; row < row_end; ++row)
        std::fill(pattern.begin() + static_cast<std::ptrdiff_t>(row * size + column_begin),
            pattern.begin() + static_cast<std::ptrdiff_t>(row * size + column_end), 1U);
}

std::pair<std::size_t, std::size_t> offset_location(
    const std::vector<std::size_t>& offsets, std::size_t index, const char* message) {
    if (offsets.empty() || index >= offsets.back()) throw std::out_of_range(message);
    const auto upper = std::upper_bound(offsets.begin(), offsets.end(), index);
    const std::size_t group = static_cast<std::size_t>(upper - offsets.begin() - 1);
    return {group, index - offsets[group]};
}

Quad4FaceBoundaryData make_boundary_data(const BoundaryConditionDefinition& boundary, bool use_displaced_geometry) {
    if (boundary.type == BoundaryConditionType::pressure)
        return {Quad4FaceBoundaryKind::pressure, CartesianTractionComponent::x, boundary.value, 0.0,
            use_displaced_geometry};
    if (boundary.type == BoundaryConditionType::traction) {
        CartesianTractionComponent component = CartesianTractionComponent::x;
        if (boundary.field == Field::displacement_y)
            component = CartesianTractionComponent::y;
        else if (boundary.field == Field::displacement_z)
            component = CartesianTractionComponent::z;
        else if (boundary.field != Field::displacement_x)
            throw std::invalid_argument("Three-dimensional traction requires a displacement field");
        return {Quad4FaceBoundaryKind::traction, component, boundary.value, 0.0, use_displaced_geometry};
    }
    return {Quad4FaceBoundaryKind::convection, CartesianTractionComponent::x, boundary.heat_transfer_coefficient,
        boundary.ambient_temperature};
}

SpatialContributionType boundary_contribution_type(BoundaryConditionType type) {
    if (type == BoundaryConditionType::pressure) return SpatialContributionType::pressure;
    if (type == BoundaryConditionType::traction) return SpatialContributionType::traction;
    return SpatialContributionType::convection;
}
} // namespace

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : SpatialLayout(definition, spatial_detail::resolve_block_ids(definition, source_mesh, true, true),
          spatial_detail::DofLayout::cartesian_3d) {
    _meshes.reserve(_block_ids.size());
    for (const std::int64_t block_id : _block_ids)
        _meshes.push_back(Hex8RegionMesh::from_unstructured_block(source_mesh, block_id));
    std::vector<std::size_t> node_counts, element_counts;
    node_counts.reserve(_meshes.size());
    element_counts.reserve(_meshes.size());
    for (const Hex8RegionMesh& mesh : _meshes) {
        node_counts.push_back(mesh.nodes().size());
        element_counts.push_back(mesh.elements().size());
    }
    initialize_counts(node_counts, element_counts);
    std::vector<std::vector<std::size_t>> region_source_node_ids;
    region_source_node_ids.reserve(_meshes.size());
    for (const Hex8RegionMesh& mesh : _meshes) region_source_node_ids.push_back(mesh.source_node_ids());
    initialize_shared_nodes(region_source_node_ids);
    _geometries.resize(_meshes.size());
    for (std::size_t region = 0; region < _meshes.size(); ++region) {
        for (const Hex8Element& element : _meshes[region].elements()) {
            Hex8Coordinates coordinates{};
            for (std::size_t node = 0; node < coordinates.size(); ++node)
                coordinates[node] = _meshes[region].nodes().at(element.nodes[node]);
            _geometries[region].push_back(make_hex8_geometry(coordinates));
        }
    }
    build_contacts(source_mesh);
    for (std::size_t boundary_index = 0; boundary_index < _definition.boundary_conditions.size(); ++boundary_index) {
        const BoundaryConditionDefinition& boundary = _definition.boundary_conditions[boundary_index];
        if (boundary.name.empty() || boundary.boundary.empty())
            throw std::invalid_argument("Cartesian three-dimensional boundary names must be nonempty");
        if (boundary.scale_with_load && !boundary.function.empty())
            throw std::invalid_argument(
                "Boundary condition cannot combine scale_with_load and a time function: " + boundary.name);
        const std::int64_t block_id = source_mesh.side_set_block_id(boundary.boundary);
        const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
        if (found == _block_ids.end())
            throw std::invalid_argument(
                "Boundary belongs to an undeclared three-dimensional block: " + boundary.boundary);
        const std::size_t region = static_cast<std::size_t>(found - _block_ids.begin());
        const Hex8RegionBoundary mapped = _meshes[region].map_side_set(source_mesh, boundary.boundary);
        if (boundary.type == BoundaryConditionType::dirichlet) {
            for (std::size_t local_node : mapped.nodes) {
                const std::size_t dof = this->dof(boundary.field, global_node(region, local_node));
                add_dirichlet(dof, boundary_index);
            }
            continue;
        }
        if (boundary.type == BoundaryConditionType::pressure || boundary.type == BoundaryConditionType::traction)
            record_configuration_warning(boundary, this->region(region));
        const std::size_t kernel = _boundary_data.size();
        const bool displaced_geometry = boundary_uses_displaced_geometry(boundary, this->region(region));
        _boundary_data.push_back(make_boundary_data(boundary, displaced_geometry));
        _boundary_definition_indices.push_back(boundary_index);
        for (const Quad4FaceElement& face : mapped.faces) {
            std::array<std::size_t, 4> nodes{};
            Quad4FaceCoordinates coordinates{};
            for (std::size_t node = 0; node < 4; ++node) nodes[node] = global_node(region, face.nodes[node]);
            for (std::size_t node = 0; node < 4; ++node)
                coordinates[node] = _meshes[region].nodes().at(face.nodes[node]);
            _boundary_contributions.push_back(
                {boundary_contribution_type(boundary.type), kernel, nodes, make_quad4_face_geometry(coordinates)});
        }
    }
    spatial_detail::validate_dirichlet_conditions(_dirichlet_conditions,
        "Cartesian three-dimensional boundary has conflicting Dirichlet values",
        "Cartesian three-dimensional boundary has duplicate Dirichlet values");
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index) {
        _kernel_data.push_back({IsotropicThermoelasticMaterial(region(region_index).material),
            region_heat_source(region_index), 0.0, region(region_index).strain_formulation});
    }
    _committed_contact_solution = initial_state();
    validate_local_state(0, contribution_count(), _committed_contact_solution);
    refresh_controls();
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh)
    : SpatialLayout(definition, spatial_detail::resolve_block_ids(definition, source_mesh, true, true),
          spatial_detail::DofLayout::cartesian_3d),
      _uses_hex20(true) {
    _hex20_meshes.reserve(_block_ids.size());
    for (const std::int64_t block_id : _block_ids)
        _hex20_meshes.push_back(Hex20RegionMesh::from_unstructured_block(source_mesh, block_id));
    std::vector<std::size_t> node_counts, element_counts;
    std::vector<std::vector<std::size_t>> region_source_node_ids;
    std::vector<std::vector<bool>> region_temperature_nodes;
    for (const Hex20RegionMesh& mesh : _hex20_meshes) {
        node_counts.push_back(mesh.nodes().size());
        element_counts.push_back(mesh.elements().size());
        region_source_node_ids.push_back(mesh.source_node_ids());
        region_temperature_nodes.push_back(mesh.temperature_nodes());
    }
    initialize_counts(node_counts, element_counts);
    initialize_mixed_shared_nodes(region_source_node_ids, region_temperature_nodes);
    _hex20_geometries.resize(_hex20_meshes.size());
    for (std::size_t region = 0; region < _hex20_meshes.size(); ++region)
        for (const Hex20Element& element : _hex20_meshes[region].elements()) {
            Hex20Coordinates coordinates{};
            for (std::size_t node = 0; node < coordinates.size(); ++node)
                coordinates[node] = _hex20_meshes[region].nodes().at(element.nodes[node]);
            _hex20_geometries[region].push_back(make_hex20_geometry(coordinates));
        }
    build_hex20_contacts(source_mesh);
    for (std::size_t boundary_index = 0; boundary_index < _definition.boundary_conditions.size(); ++boundary_index) {
        const BoundaryConditionDefinition& boundary = _definition.boundary_conditions[boundary_index];
        if (boundary.name.empty() || boundary.boundary.empty())
            throw std::invalid_argument("HEX20 boundary names must be nonempty");
        if (boundary.scale_with_load && !boundary.function.empty())
            throw std::invalid_argument(
                "Boundary condition cannot combine scale_with_load and a time function: " + boundary.name);
        const std::int64_t block_id = source_mesh.side_set_block_id(boundary.boundary);
        const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
        if (found == _block_ids.end())
            throw std::invalid_argument("Boundary belongs to an undeclared HEX20 block: " + boundary.boundary);
        const std::size_t region = static_cast<std::size_t>(found - _block_ids.begin());
        const Hex20RegionBoundary mapped = _hex20_meshes[region].map_side_set(source_mesh, boundary.boundary);
        if (boundary.type == BoundaryConditionType::dirichlet) {
            const std::vector<std::size_t>& nodes =
                boundary.field == Field::temperature ? mapped.temperature_nodes : mapped.displacement_nodes;
            for (const std::size_t local_node : nodes) {
                const std::size_t global = boundary.field == Field::temperature
                                               ? global_temperature_node(region, local_node)
                                               : global_node(region, local_node);
                add_dirichlet(dof(boundary.field, global), boundary_index);
            }
            continue;
        }
        if (boundary.type == BoundaryConditionType::pressure || boundary.type == BoundaryConditionType::traction)
            record_configuration_warning(boundary, this->region(region));
        const std::size_t kernel = _boundary_data.size();
        const bool displaced_geometry = boundary_uses_displaced_geometry(boundary, this->region(region));
        _boundary_data.push_back(make_boundary_data(boundary, displaced_geometry));
        _boundary_definition_indices.push_back(boundary_index);
        for (const Quad8FaceElement& face : mapped.faces) {
            std::array<std::size_t, 4> temperature_nodes{};
            std::array<std::size_t, 8> displacement_nodes{};
            Quad8FaceCoordinates coordinates{};
            for (std::size_t node = 0; node < 8; ++node) {
                displacement_nodes[node] = global_node(region, face.nodes[node]);
                coordinates[node] = _hex20_meshes[region].nodes().at(face.nodes[node]);
                if (node < 4) temperature_nodes[node] = global_temperature_node(region, face.nodes[node]);
            }
            _hex20_boundary_contributions.push_back({boundary_contribution_type(boundary.type), kernel,
                temperature_nodes, displacement_nodes, make_quad8_face_geometry(coordinates)});
        }
    }
    spatial_detail::validate_dirichlet_conditions(_dirichlet_conditions,
        "HEX20 boundary has conflicting Dirichlet values", "HEX20 boundary has duplicate Dirichlet values");
    for (std::size_t region = 0; region < region_count(); ++region)
        _kernel_data.push_back({IsotropicThermoelasticMaterial(this->region(region).material),
            region_heat_source(region), 0.0, this->region(region).strain_formulation});
    _committed_contact_solution = initial_state();
    validate_local_state(0, contribution_count(), _committed_contact_solution);
    refresh_controls();
}

SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    const ContributionRanges ranges = contribution_ranges();
    if (index >= ranges.end) throw std::out_of_range("Three-dimensional contribution index is out of range");
    if (index < ranges.thermal_begin) return SpatialContributionType::volume;
    if (index < ranges.mechanical_begin) return SpatialContributionType::thermal_contact;
    if (index < ranges.boundary_begin) return SpatialContributionType::mechanical_contact;
    return _uses_hex20 ? _hex20_boundary_contributions.at(index - ranges.boundary_begin).type
                       : _boundary_contributions.at(index - ranges.boundary_begin).type;
}

const Hex8Geometry& SpatialAssembly::region_element_geometry(std::size_t region, std::size_t element_index) const {
    return _geometries.at(region).at(element_index);
}

const Hex20Geometry& SpatialAssembly::hex20_region_element_geometry(
    std::size_t region, std::size_t element_index) const {
    return _hex20_geometries.at(region).at(element_index);
}

void SpatialAssembly::set_load_factor(double value) {
    set_load_factor_value(value);
    refresh_controls();
}

void SpatialAssembly::set_time(double value) {
    set_time_value(value);
    for (CartesianThermoelasticData& kernel_data : _kernel_data) kernel_data.time = value;
    refresh_controls();
}

void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional state size does not match the problem");
    if (!std::all_of(state.begin(), state.end(), [](double value) { return std::isfinite(value); }))
        throw std::domain_error("Three-dimensional state must contain only finite values");
    for (std::size_t region_index = 0; region_index < region_count(); ++region_index) {
        if (region(region_index).strain_formulation != StrainFormulation::finite) continue;
        for (std::size_t element = 0; element < region_element_count(region_index); ++element) {
            if (_uses_hex20) {
                const Hex20LocalValues local = hex20_volume_state(region_element_offset(region_index) + element, state);
                for (const Hex20MechanicalQuadraturePoint& point :
                    hex20_region_element_geometry(region_index, element).mechanical_points)
                    validate_hex20_deformation(point, local);
                continue;
            }
            const Hex8LocalValues local = volume_state(region_element_offset(region_index) + element, state);
            for (const Hex8QuadraturePoint& point : region_element_geometry(region_index, element).points)
                validate_cartesian_deformation(point, local);
        }
    }
    validate_local_state(0, contribution_count(), state);
}

void SpatialAssembly::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    if (index < volume_contribution_count()) {
        const auto location = element_location(index);
        if (_uses_hex20) {
            const Hex20Element& element = _hex20_meshes[location.first].elements()[location.second];
            Hex20LocalDofs fixed{};
            for (std::size_t node = 0; node < 8; ++node)
                fixed[node] = dof(Field::temperature, global_temperature_node(location.first, element.nodes[node]));
            for (std::size_t node = 0; node < 20; ++node) {
                const std::size_t global = global_node(location.first, element.nodes[node]);
                fixed[8 + node] = dof(Field::displacement_x, global);
                fixed[28 + node] = dof(Field::displacement_y, global);
                fixed[48 + node] = dof(Field::displacement_z, global);
            }
            dofs.assign(fixed.begin(), fixed.end());
            return;
        }
        std::array<std::size_t, 8> nodes{};
        const Hex8Element& element = _meshes[location.first].elements()[location.second];
        for (std::size_t node = 0; node < 8; ++node) nodes[node] = global_node(location.first, element.nodes[node]);
        Hex8LocalDofs fixed{};
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            fixed[node] = dof(Field::temperature, nodes[node]);
            fixed[8 + node] = dof(Field::displacement_x, nodes[node]);
            fixed[16 + node] = dof(Field::displacement_y, nodes[node]);
            fixed[24 + node] = dof(Field::displacement_z, nodes[node]);
        }
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    const ContributionRanges ranges = contribution_ranges();
    if (index < ranges.mechanical_begin) {
        const std::size_t point = index - ranges.thermal_begin;
        if (_uses_hex20) {
            const Quad8SurfaceContactLocalDofs fixed =
                hex20_contact_dofs(hex20_thermal_candidate(point, _thermal_active_primary.at(point)));
            dofs.assign(fixed.begin(), fixed.end());
            return;
        }
        const Quad4SurfaceContactLocalDofs fixed =
            contact_dofs(thermal_candidate(point, _thermal_active_primary.at(point)).nodes);
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    if (index < ranges.boundary_begin) {
        const std::size_t point = index - ranges.mechanical_begin;
        if (_uses_hex20) {
            if (point >= _hex20_mechanical_points.size()) {
                hex20_averaged_constraint_dofs(
                    _hex20_averaged_constraints.at(point - _hex20_mechanical_points.size()), dofs);
                return;
            }
            const Quad8SurfaceContactLocalDofs fixed =
                hex20_contact_dofs(hex20_mechanical_candidate(point, _mechanical_active_primary.at(point)));
            dofs.assign(fixed.begin(), fixed.end());
            return;
        }
        const Quad4SurfaceContactLocalDofs fixed =
            contact_dofs(mechanical_candidate(point, _mechanical_active_primary.at(point)).nodes);
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    if (_uses_hex20) {
        const Hex20BoundaryContribution& entry = _hex20_boundary_contributions.at(index - ranges.boundary_begin);
        Quad8FaceLocalDofs fixed{};
        for (std::size_t node = 0; node < 4; ++node)
            fixed[node] = dof(Field::temperature, entry.temperature_nodes[node]);
        for (std::size_t node = 0; node < 8; ++node) {
            fixed[4 + node] = dof(Field::displacement_x, entry.displacement_nodes[node]);
            fixed[12 + node] = dof(Field::displacement_y, entry.displacement_nodes[node]);
            fixed[20 + node] = dof(Field::displacement_z, entry.displacement_nodes[node]);
        }
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    const BoundaryContribution& entry = _boundary_contributions.at(index - ranges.boundary_begin);
    Quad4FaceLocalDofs fixed{};
    for (std::size_t node = 0; node < entry.nodes.size(); ++node) {
        fixed[node] = dof(Field::temperature, entry.nodes[node]);
        fixed[4 + node] = dof(Field::displacement_x, entry.nodes[node]);
        fixed[8 + node] = dof(Field::displacement_y, entry.nodes[node]);
        fixed[12 + node] = dof(Field::displacement_z, entry.nodes[node]);
    }
    dofs.assign(fixed.begin(), fixed.end());
}

void SpatialAssembly::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    const ContributionRanges ranges = contribution_ranges();
    if (index >= ranges.end) throw std::out_of_range("Three-dimensional contribution index is out of range");
    if (index < ranges.thermal_begin) {
        if (_uses_hex20) {
            pattern.assign(hex20_local_dof_count * hex20_local_dof_count, 0U);
            set_pattern_block(pattern, hex20_local_dof_count, 0, 8, 0, 8);
            set_pattern_block(pattern, hex20_local_dof_count, 8, hex20_local_dof_count, 0, hex20_local_dof_count);
            return;
        }
        pattern.assign(hex8_local_dof_count * hex8_local_dof_count, 0U);
        set_pattern_block(pattern, hex8_local_dof_count, 0, 8, 0, 8);
        set_pattern_block(pattern, hex8_local_dof_count, 8, hex8_local_dof_count, 0, hex8_local_dof_count);
        return;
    }
    if (index < ranges.mechanical_begin) {
        if (_uses_hex20) {
            pattern.assign(quad8_surface_contact_local_dof_count * quad8_surface_contact_local_dof_count, 0U);
            set_pattern_block(
                pattern, quad8_surface_contact_local_dof_count, 0, 8, 0, quad8_surface_contact_local_dof_count);
            return;
        }
        pattern.assign(quad4_surface_contact_local_dof_count * quad4_surface_contact_local_dof_count, 0U);
        set_pattern_block(
            pattern, quad4_surface_contact_local_dof_count, 0, 8, 0, quad4_surface_contact_local_dof_count);
        return;
    }
    if (index < ranges.boundary_begin) {
        if (_uses_hex20) {
            const std::size_t point = index - ranges.mechanical_begin;
            if (point >= _hex20_mechanical_points.size()) {
                const std::size_t local_size =
                    3 * _hex20_averaged_constraints.at(point - _hex20_mechanical_points.size()).nodes.size();
                pattern.assign(local_size * local_size, 1U);
                return;
            }
            pattern.assign(quad8_surface_contact_local_dof_count * quad8_surface_contact_local_dof_count, 0U);
            set_pattern_block(pattern, quad8_surface_contact_local_dof_count, 8, quad8_surface_contact_local_dof_count,
                8, quad8_surface_contact_local_dof_count);
            return;
        }
        pattern.assign(quad4_surface_contact_local_dof_count * quad4_surface_contact_local_dof_count, 0U);
        set_pattern_block(pattern, quad4_surface_contact_local_dof_count, 8, quad4_surface_contact_local_dof_count, 8,
            quad4_surface_contact_local_dof_count);
        return;
    }
    const std::size_t kernel = _uses_hex20 ? _hex20_boundary_contributions.at(index - ranges.boundary_begin).kernel
                                           : _boundary_contributions.at(index - ranges.boundary_begin).kernel;
    const Quad4FaceBoundaryData& data = _boundary_data.at(kernel);
    if (_uses_hex20) {
        pattern.assign(quad8_face_local_dof_count * quad8_face_local_dof_count, 0U);
        if (data.kind == Quad4FaceBoundaryKind::convection) {
            set_pattern_block(pattern, quad8_face_local_dof_count, 0, 4, 0, 4);
        } else if (data.use_displaced_geometry) {
            const std::size_t row_begin = data.kind == Quad4FaceBoundaryKind::pressure
                                              ? 4
                                              : (data.component == CartesianTractionComponent::x
                                                        ? 4
                                                        : (data.component == CartesianTractionComponent::y ? 12 : 20));
            const std::size_t row_end = data.kind == Quad4FaceBoundaryKind::pressure ? 28 : row_begin + 8;
            set_pattern_block(pattern, quad8_face_local_dof_count, row_begin, row_end, 4, 28);
        }
        return;
    }
    pattern.assign(quad4_face_local_dof_count * quad4_face_local_dof_count, 0U);
    if (data.kind == Quad4FaceBoundaryKind::convection) {
        set_pattern_block(pattern, quad4_face_local_dof_count, 0, 4, 0, 4);
    } else if (data.use_displaced_geometry) {
        const std::size_t row_begin = data.kind == Quad4FaceBoundaryKind::pressure
                                          ? 4
                                          : (data.component == CartesianTractionComponent::x
                                                    ? 4
                                                    : (data.component == CartesianTractionComponent::y ? 8 : 12));
        const std::size_t row_end = data.kind == Quad4FaceBoundaryKind::pressure ? 16 : row_begin + 4;
        set_pattern_block(pattern, quad4_face_local_dof_count, row_begin, row_end, 4, 16);
    }
}

void SpatialAssembly::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    if (index < volume_contribution_count()) {
        contribution_dofs(index, dofs);
        return;
    }
    index -= volume_contribution_count();
    if (_uses_hex20) {
        if (index >= _sparsity_contact_offsets.back()) {
            hex20_averaged_constraint_dofs(
                _hex20_averaged_constraints.at(index - _sparsity_contact_offsets.back()), dofs);
            return;
        }
        const Hex20SparsityContact contact = hex20_sparsity_contact(index);
        Quad8SurfaceContactLocalDofs fixed{};
        for (std::size_t node = 0; node < 4; ++node) {
            fixed[node] = dof(Field::temperature, contact.temperature_nodes[node]);
            fixed[4 + node] = dof(Field::temperature, contact.primary_temperature_nodes[node]);
        }
        for (std::size_t component = 0; component < 3; ++component) {
            const Field field = component == 0 ? Field::displacement_x
                                               : (component == 1 ? Field::displacement_y : Field::displacement_z);
            const std::size_t offset = 8 + 16 * component;
            for (std::size_t node = 0; node < 8; ++node) {
                fixed[offset + node] = dof(field, contact.displacement_nodes[node]);
                fixed[offset + 8 + node] = dof(field, contact.primary_displacement_nodes[node]);
            }
        }
        dofs.assign(fixed.begin(), fixed.end());
        return;
    }
    const Quad4SurfaceContactLocalDofs fixed = contact_dofs(sparsity_contact(index).nodes);
    dofs.assign(fixed.begin(), fixed.end());
}

void SpatialAssembly::sparsity_contribution_jacobian_pattern(
    std::size_t index, std::vector<unsigned char>& pattern) const {
    if (index < volume_contribution_count()) {
        contribution_jacobian_pattern(index, pattern);
        return;
    }
    index -= volume_contribution_count();
    if (_uses_hex20) {
        if (index >= _sparsity_contact_offsets.back()) {
            const std::size_t local_size =
                3 * _hex20_averaged_constraints.at(index - _sparsity_contact_offsets.back()).nodes.size();
            pattern.assign(local_size * local_size, 1U);
            return;
        }
        const Hex20SparsityContact contact = hex20_sparsity_contact(index);
        pattern.assign(quad8_surface_contact_local_dof_count * quad8_surface_contact_local_dof_count, 0U);
        if (contact.thermal)
            set_pattern_block(
                pattern, quad8_surface_contact_local_dof_count, 0, 8, 0, quad8_surface_contact_local_dof_count);
        if (contact.mechanical)
            set_pattern_block(pattern, quad8_surface_contact_local_dof_count, 8, quad8_surface_contact_local_dof_count,
                8, quad8_surface_contact_local_dof_count);
        return;
    }
    const SparsityContact contact = sparsity_contact(index);
    pattern.assign(quad4_surface_contact_local_dof_count * quad4_surface_contact_local_dof_count, 0U);
    if (contact.thermal)
        set_pattern_block(
            pattern, quad4_surface_contact_local_dof_count, 0, 8, 0, quad4_surface_contact_local_dof_count);
    if (contact.mechanical)
        set_pattern_block(pattern, quad4_surface_contact_local_dof_count, 8, quad4_surface_contact_local_dof_count, 8,
            quad4_surface_contact_local_dof_count);
}

Hex8LocalValues SpatialAssembly::volume_state(std::size_t index, const std::vector<double>& global_state) const {
    if (_uses_hex20) throw std::logic_error("HEX8 volume state requested from a HEX20 problem");
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    if (dofs.size() != hex8_local_dof_count)
        throw std::logic_error("HEX8 volume contribution has an invalid DOF layout");
    Hex8LocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

Hex20LocalValues SpatialAssembly::hex20_volume_state(std::size_t index, const std::vector<double>& global_state) const {
    if (!_uses_hex20) throw std::logic_error("HEX20 volume state requested from a HEX8 problem");
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    if (dofs.size() != hex20_local_dof_count)
        throw std::logic_error("HEX20 volume contribution has an invalid DOF layout");
    Hex20LocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

void SpatialAssembly::compute_contribution(std::size_t index, const std::vector<double>& state,
    const std::vector<double>* committed_solution, const CartesianMaterialHistory* committed_material, double time_step,
    std::vector<double>& residual, std::vector<double>* jacobian, bool include_thermal_time_term) const {
    if (index < volume_contribution_count()) {
        if (_uses_hex20) {
            if (state.size() != hex20_local_dof_count)
                throw std::invalid_argument("HEX20 contribution state must contain 68 DOFs");
            const auto location = element_location(index);
            Hex20LocalValues current{};
            std::copy(state.begin(), state.end(), current.begin());
            const Hex20LocalValues committed =
                committed_solution == nullptr ? Hex20LocalValues{} : hex20_volume_state(index, *committed_solution);
            Hex20LocalJacobian local_jacobian{};
            const Hex20LocalResidual result =
                committed_material == nullptr
                    ? compute_hex20_thermoelastic(_kernel_data[location.first],
                          hex20_region_element_geometry(location.first, location.second), current,
                          committed_solution == nullptr ? nullptr : &committed, time_step,
                          jacobian == nullptr ? nullptr : &local_jacobian)
                    : compute_hex20_transient(_kernel_data[location.first],
                          hex20_region_element_geometry(location.first, location.second), current, committed,
                          *committed_material, time_step, jacobian == nullptr ? nullptr : &local_jacobian,
                          include_thermal_time_term);
            residual.assign(result.begin(), result.end());
            if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
            return;
        }
        if (state.size() != hex8_local_dof_count)
            throw std::invalid_argument("HEX8 contribution state must contain 32 DOFs");
        const auto location = element_location(index);
        Hex8LocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        const Hex8LocalValues committed =
            committed_solution == nullptr ? Hex8LocalValues{} : volume_state(index, *committed_solution);
        Hex8LocalJacobian local_jacobian{};
        const Hex8LocalResidual result =
            committed_material == nullptr
                ? compute_hex8_thermoelastic(_kernel_data[location.first],
                      region_element_geometry(location.first, location.second), current,
                      committed_solution == nullptr ? nullptr : &committed, time_step,
                      jacobian == nullptr ? nullptr : &local_jacobian)
                : compute_hex8_transient(_kernel_data[location.first],
                      region_element_geometry(location.first, location.second), current, committed, *committed_material,
                      time_step, jacobian == nullptr ? nullptr : &local_jacobian, include_thermal_time_term);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    const ContributionRanges ranges = contribution_ranges();
    if (index < ranges.mechanical_begin) {
        if (_uses_hex20) {
            if (state.size() != quad8_surface_contact_local_dof_count)
                throw std::invalid_argument("HEX20 thermal-contact state must contain 56 DOFs");
            const std::size_t point = index - ranges.thermal_begin;
            const Hex20ThermalCandidate entry = hex20_thermal_candidate(point, _thermal_active_primary.at(point));
            Quad8SurfaceContactLocalValues current{};
            std::copy(state.begin(), state.end(), current.begin());
            Quad8SurfaceContactLocalJacobian local_jacobian{};
            const Quad8SurfaceContactLocalResidual result =
                compute_quad8_to_quad8_gap_heat(_thermal_properties[entry.contact], entry.geometry, current,
                    jacobian == nullptr ? nullptr : &local_jacobian);
            residual.assign(result.begin(), result.end());
            if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
            return;
        }
        if (state.size() != quad4_surface_contact_local_dof_count)
            throw std::invalid_argument("Three-dimensional thermal-contact state must contain 32 DOFs");
        const std::size_t point = index - ranges.thermal_begin;
        const ThermalCandidate entry = thermal_candidate(point, _thermal_active_primary.at(point));
        Quad4SurfaceContactLocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        Quad4SurfaceContactLocalJacobian local_jacobian{};
        const Quad4SurfaceContactLocalResidual result =
            compute_quad4_to_quad4_gap_heat(_thermal_properties[entry.contact], entry.geometry, current,
                jacobian == nullptr ? nullptr : &local_jacobian);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    if (index < ranges.boundary_begin) {
        if (_uses_hex20) {
            const std::size_t point = index - ranges.mechanical_begin;
            if (point >= _hex20_mechanical_points.size()) {
                compute_hex20_averaged_constraint(
                    _hex20_averaged_constraints.at(point - _hex20_mechanical_points.size()), state, residual, jacobian);
                return;
            }
            if (state.size() != quad8_surface_contact_local_dof_count)
                throw std::invalid_argument("HEX20 mechanical-contact state must contain 56 DOFs");
            const Hex20MechanicalCandidate entry =
                hex20_mechanical_candidate(point, _mechanical_active_primary.at(point));
            Quad8SurfaceContactLocalValues current{};
            std::copy(state.begin(), state.end(), current.begin());
            const Quad8SurfaceContactLocalValues committed = hex20_contact_state(entry, _committed_contact_solution);
            Quad8SurfaceContactLocalJacobian local_jacobian{};
            const Quad8SurfaceContactLocalResidual result =
                entry.surface_to_surface
                    ? compute_quad8_to_quad8_contact(_mechanical_properties[entry.contact], entry.surface_geometry,
                          current, committed, _contact_histories[entry.contact][entry.secondary],
                          jacobian == nullptr ? nullptr : &local_jacobian)
                    : compute_node_to_quad8_contact(_mechanical_properties[entry.contact], entry.node_geometry, current,
                          committed, _contact_histories[entry.contact][entry.secondary],
                          jacobian == nullptr ? nullptr : &local_jacobian);
            residual.assign(result.begin(), result.end());
            if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
            return;
        }
        if (state.size() != quad4_surface_contact_local_dof_count)
            throw std::invalid_argument("Three-dimensional mechanical-contact state must contain 32 DOFs");
        const std::size_t point = index - ranges.mechanical_begin;
        const MechanicalCandidate entry = mechanical_candidate(point, _mechanical_active_primary.at(point));
        Quad4SurfaceContactLocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        const Quad4SurfaceContactLocalValues committed = contribution_state(index, _committed_contact_solution);
        Quad4SurfaceContactLocalJacobian local_jacobian{};
        const Quad4SurfaceContactLocalResidual result =
            compute_node_to_quad4_contact(_mechanical_properties[entry.contact], entry.geometry, current, committed,
                _contact_histories[entry.contact][entry.secondary], jacobian == nullptr ? nullptr : &local_jacobian);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    if (_uses_hex20) {
        if (state.size() != quad8_face_local_dof_count)
            throw std::invalid_argument("Three-dimensional quadratic face state must contain 28 DOFs");
        const Hex20BoundaryContribution& entry = _hex20_boundary_contributions.at(index - ranges.boundary_begin);
        Quad8FaceLocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        Quad8FaceLocalJacobian local_jacobian{};
        const Quad8FaceLocalResidual result = compute_quad8_face_boundary(
            _boundary_data[entry.kernel], entry.geometry, current, jacobian == nullptr ? nullptr : &local_jacobian);
        residual.assign(result.begin(), result.end());
        if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
        return;
    }
    if (state.size() != quad4_face_local_dof_count)
        throw std::invalid_argument("Three-dimensional face state must contain 16 DOFs");
    const BoundaryContribution& entry = _boundary_contributions.at(index - ranges.boundary_begin);
    Quad4FaceLocalValues current{};
    std::copy(state.begin(), state.end(), current.begin());
    Quad4FaceLocalJacobian local_jacobian{};
    const Quad4FaceLocalResidual result = compute_quad4_face_boundary(
        _boundary_data[entry.kernel], entry.geometry, current, jacobian == nullptr ? nullptr : &local_jacobian);
    residual.assign(result.begin(), result.end());
    if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
}

CartesianMaterialHistory SpatialAssembly::transient_update(std::size_t region, std::size_t element,
    const Hex8LocalValues& state, const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material, double time_step) const {
    return compute_hex8_transient_update(_kernel_data.at(region), region_element_geometry(region, element), state,
        committed_state, committed_material, time_step);
}

CartesianMaterialHistory SpatialAssembly::transient_update(std::size_t region, std::size_t element,
    const Hex20LocalValues& state, const Hex20LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material, double time_step) const {
    return compute_hex20_transient_update(_kernel_data.at(region), hex20_region_element_geometry(region, element),
        state, committed_state, committed_material, time_step);
}

std::array<SymmetricTensor3Values, 8> SpatialAssembly::stress(
    std::size_t region, std::size_t element, const std::vector<double>& state) const {
    return compute_hex8_stress(_kernel_data.at(region), region_element_geometry(region, element),
        volume_state(region_element_offset(region) + element, state));
}

std::array<SymmetricTensor3Values, 27> SpatialAssembly::hex20_stress(
    std::size_t region, std::size_t element, const std::vector<double>& state) const {
    return compute_hex20_stress(_kernel_data.at(region), hex20_region_element_geometry(region, element),
        hex20_volume_state(region_element_offset(region) + element, state));
}

double SpatialAssembly::heat_capacity(std::size_t region, double temperature, const CartesianPoint3& position) const {
    const CartesianThermoelasticData& data = _kernel_data.at(region);
    return data.material.heat_capacity(temperature, {data.time, position.x, position.y, position.z}).value();
}

SpatialAssembly::ContributionRanges SpatialAssembly::contribution_ranges() const noexcept {
    const std::size_t thermal_begin = volume_contribution_count(),
                      mechanical_begin = thermal_begin + _thermal_contact_offsets.back(),
                      mechanical_count = _uses_hex20
                                             ? _hex20_mechanical_points.size() + _hex20_averaged_constraints.size()
                                             : _mechanical_points.size(),
                      boundary_begin = mechanical_begin + mechanical_count;
    const std::size_t boundary_count =
        _uses_hex20 ? _hex20_boundary_contributions.size() : _boundary_contributions.size();
    return {thermal_begin, mechanical_begin, boundary_begin, boundary_begin + boundary_count};
}

std::size_t SpatialAssembly::sparsity_contribution_count() const noexcept {
    // Every face boundary block is a subset of its adjacent volume block. Contact quadrature points and face nodes
    // that share one secondary-face/primary-face pair also have the same 32-DOF graph, so one representative preserves
    // the complete graph without repeating it for each runtime contribution.
    return volume_contribution_count() + _sparsity_contact_offsets.back() +
           (_uses_hex20 ? _hex20_averaged_constraints.size() : 0);
}

std::size_t SpatialAssembly::contribution_work(std::size_t index, std::size_t partition_count) const {
    // These relative units track the measured AD-local work of eight-point finite-strain volume integration and the
    // three surface kernels. They affect only the contiguous MPI ownership boundary, never the residual or Jacobian.
    const ContributionRanges ranges = contribution_ranges();
    if (index < ranges.thermal_begin) return _uses_hex20 ? 216 : 64;
    if (partition_count < 4) {
        if (index < ranges.mechanical_begin) return 10;
        if (index < ranges.boundary_begin) return 11;
        if (index < ranges.end) return 1;
    } else {
        // With four or more partitions the last contiguous interval owns every surface contribution. Account for
        // residual work as well as AD Jacobian work so that this interval does not become the synchronization tail.
        if (index < ranges.mechanical_begin) return 13;
        if (index < ranges.boundary_begin) return 14;
        if (index < ranges.end) return 2;
    }
    throw std::out_of_range("Three-dimensional contribution work index is out of range");
}

std::pair<std::size_t, std::size_t> SpatialAssembly::contribution_partition(
    std::size_t partition, std::size_t partition_count) const {
    if (partition_count == 0 || partition >= partition_count)
        throw std::out_of_range("Three-dimensional contribution partition is out of range");
    const std::size_t count = contribution_count();
    std::size_t total_work = 0;
    for (std::size_t entry = 0; entry < count; ++entry) {
        const std::size_t work = contribution_work(entry, partition_count);
        if (total_work > std::numeric_limits<std::size_t>::max() - work)
            throw std::overflow_error("Three-dimensional contribution work exceeds size_t range");
        total_work += work;
    }
    const auto boundary = [&](std::size_t boundary_partition) {
        if (boundary_partition == 0) return std::size_t{0};
        if (boundary_partition == partition_count) return count;
        const std::size_t target = (total_work / partition_count) * boundary_partition +
                                   ((total_work % partition_count) * boundary_partition) / partition_count;
        std::size_t accumulated = 0;
        for (std::size_t entry = 0; entry < count; ++entry) {
            const std::size_t next = accumulated + contribution_work(entry, partition_count);
            if (next >= target) return target - accumulated < next - target ? entry : entry + 1U;
            accumulated = next;
        }
        return count;
    };
    return {boundary(partition), boundary(partition + 1U)};
}

ResolvedBoundary SpatialAssembly::resolve_boundary(
    const UnstructuredHex8Mesh& source_mesh, const std::string& name) const {
    const std::int64_t block_id = source_mesh.side_set_block_id(name);
    const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
    if (found == _block_ids.end())
        throw std::invalid_argument("Contact boundary belongs to an undeclared three-dimensional block: " + name);
    const std::size_t region = static_cast<std::size_t>(found - _block_ids.begin());
    Hex8RegionBoundary boundary = _meshes[region].map_side_set(source_mesh, name);
    if (boundary.faces.empty()) throw std::invalid_argument("Three-dimensional contact side set is empty: " + name);
    return {region, std::move(boundary)};
}

ResolvedHex20Boundary SpatialAssembly::resolve_boundary(
    const UnstructuredHex20Mesh& source_mesh, const std::string& name) const {
    const std::int64_t block_id = source_mesh.side_set_block_id(name);
    const auto found = std::find(_block_ids.begin(), _block_ids.end(), block_id);
    if (found == _block_ids.end())
        throw std::invalid_argument("Contact boundary belongs to an undeclared HEX20 block: " + name);
    const std::size_t region = static_cast<std::size_t>(found - _block_ids.begin());
    Hex20RegionBoundary boundary = _hex20_meshes[region].map_side_set(source_mesh, name);
    if (boundary.faces.empty()) throw std::invalid_argument("HEX20 contact side set is empty: " + name);
    return {region, std::move(boundary)};
}

Quad4SurfaceContactLocalDofs SpatialAssembly::contact_dofs(const std::array<std::size_t, 8>& nodes) const {
    Quad4SurfaceContactLocalDofs result{};
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        result[node] = dof(Field::temperature, nodes[node]);
        result[8 + node] = dof(Field::displacement_x, nodes[node]);
        result[16 + node] = dof(Field::displacement_y, nodes[node]);
        result[24 + node] = dof(Field::displacement_z, nodes[node]);
    }
    return result;
}

Quad4SurfaceContactLocalValues SpatialAssembly::contribution_state(
    std::size_t index, const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional contact global state has the wrong size");
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    if (dofs.size() != quad4_surface_contact_local_dof_count)
        throw std::logic_error("Three-dimensional contact contribution has an invalid DOF layout");
    Quad4SurfaceContactLocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

Quad4SurfaceContactLocalValues SpatialAssembly::contact_state(
    const std::array<std::size_t, 8>& nodes, const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional contact state has the wrong size");
    const Quad4SurfaceContactLocalDofs dofs = contact_dofs(nodes);
    Quad4SurfaceContactLocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

Quad8SurfaceContactLocalDofs SpatialAssembly::hex20_contact_dofs(const Hex20ThermalCandidate& candidate) const {
    Quad8SurfaceContactLocalDofs result{};
    for (std::size_t node = 0; node < 4; ++node) {
        result[node] = dof(Field::temperature, candidate.secondary_temperature_nodes[node]);
        result[4 + node] = dof(Field::temperature, candidate.primary_temperature_nodes[node]);
    }
    const std::array<Field, 3> displacement_fields = {
        Field::displacement_x, Field::displacement_y, Field::displacement_z};
    for (std::size_t component = 0; component < 3; ++component) {
        const std::size_t offset = 8 + 16 * component;
        for (std::size_t node = 0; node < 8; ++node) {
            result[offset + node] = dof(displacement_fields[component], candidate.secondary_displacement_nodes[node]);
            result[offset + 8 + node] = dof(displacement_fields[component], candidate.primary_displacement_nodes[node]);
        }
    }
    return result;
}

Quad8SurfaceContactLocalDofs SpatialAssembly::hex20_contact_dofs(const Hex20MechanicalCandidate& candidate) const {
    Hex20ThermalCandidate thermal{candidate.contact, candidate.secondary_temperature_nodes,
        candidate.primary_temperature_nodes, candidate.secondary_displacement_nodes,
        candidate.primary_displacement_nodes, {}, candidate.primary};
    return hex20_contact_dofs(thermal);
}

void SpatialAssembly::hex20_averaged_constraint_dofs(
    const Hex20AveragedConstraint& constraint, std::vector<std::size_t>& dofs) const {
    dofs.clear();
    dofs.reserve(3 * constraint.nodes.size());
    for (const Field field : {Field::displacement_x, Field::displacement_y, Field::displacement_z})
        for (const std::size_t node : constraint.nodes) dofs.push_back(dof(field, node));
}

SpatialAssembly::Hex20AveragedConstraintValue SpatialAssembly::hex20_averaged_constraint_value(
    const Hex20AveragedConstraint& constraint, const std::vector<double>& state,
    const std::vector<double>& committed_state, const ContactPointHistory& history) const {
    const std::size_t node_count = constraint.nodes.size();
    if (state.size() != 3 * node_count || committed_state.size() != state.size())
        throw std::invalid_argument("HEX20 averaged-contact state has the wrong size");
    std::array<double, 3> relative{}, committed_relative{};
    for (std::size_t component = 0; component < 3; ++component) {
        for (std::size_t node = 0; node < node_count; ++node) {
            const double coefficient = constraint.gap_coefficients[node];
            relative[component] += coefficient * state[component * node_count + node];
            committed_relative[component] += coefficient * committed_state[component * node_count + node];
        }
    }
    const std::array<double, 3> normal = {constraint.normal.x, constraint.normal.y, constraint.normal.z};
    Hex20AveragedConstraintValue result{};
    result.gap = constraint.reference_gap;
    for (std::size_t component = 0; component < 3; ++component) result.gap += normal[component] * relative[component];
    double normal_relative = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        normal_relative += normal[component] * relative[component];
    for (std::size_t component = 0; component < 3; ++component)
        result.tangential_slip[component] = -relative[component] + normal_relative * normal[component];
    const NormalContactProperties& properties = _mechanical_properties[constraint.contact];
    result.pressure = std::max(-properties.penalty * result.gap, 0.0);
    result.force = result.pressure * constraint.area;
    if (properties.friction_coefficient == 0.0 || !(result.pressure > 0.0)) return result;

    std::array<double, 3> relative_increment{};
    double normal_increment = 0.0;
    for (std::size_t component = 0; component < 3; ++component) {
        // The normal gap is primary minus secondary, whereas the established
        // friction history convention is secondary minus primary.
        relative_increment[component] = committed_relative[component] - relative[component];
        normal_increment += normal[component] * relative_increment[component];
    }
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] = history.cartesian_elastic_tangential_slip[component] +
                                                    relative_increment[component] -
                                                    normal_increment * normal[component];
    double history_normal = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        history_normal += result.elastic_tangential_slip[component] * normal[component];
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] -= history_normal * normal[component];
    const double sliding_limit = properties.friction_coefficient * result.pressure;
    result.stick_stiffness =
        properties.maximum_elastic_slip > 0.0 ? sliding_limit / properties.maximum_elastic_slip : properties.penalty;
    if (!std::isfinite(result.stick_stiffness) || !(result.stick_stiffness > 0.0))
        throw std::domain_error("HEX20 averaged contact has a nonpositive tangential stick stiffness");
    for (std::size_t component = 0; component < 3; ++component) {
        result.trial_tangential_traction[component] =
            result.stick_stiffness * result.elastic_tangential_slip[component];
    }
    result.trial_tangential_magnitude = std::hypot(
        result.trial_tangential_traction[0], result.trial_tangential_traction[1], result.trial_tangential_traction[2]);
    if (result.trial_tangential_magnitude < sliding_limit ||
        (result.trial_tangential_magnitude == sliding_limit && !history.sliding)) {
        result.tangential_traction = result.trial_tangential_traction;
    } else {
        if (!(result.trial_tangential_magnitude > 0.0))
            throw std::domain_error("HEX20 averaged sliding contact has an undefined tangential direction");
        for (std::size_t component = 0; component < 3; ++component) {
            result.tangential_traction[component] =
                sliding_limit * result.trial_tangential_traction[component] / result.trial_tangential_magnitude;
            result.elastic_tangential_slip[component] = result.tangential_traction[component] / result.stick_stiffness;
        }
        result.sliding = true;
    }
    result.tangential_force = constraint.area * std::hypot(result.tangential_traction[0], result.tangential_traction[1],
                                                    result.tangential_traction[2]);
    return result;
}

void SpatialAssembly::compute_hex20_averaged_constraint(const Hex20AveragedConstraint& constraint,
    const std::vector<double>& state, std::vector<double>& residual, std::vector<double>* jacobian) const {
    const std::size_t node_count = constraint.nodes.size(), local_size = 3 * node_count;
    std::vector<std::size_t> dofs;
    hex20_averaged_constraint_dofs(constraint, dofs);
    std::vector<double> committed_state(local_size);
    for (std::size_t local = 0; local < local_size; ++local)
        committed_state[local] = _committed_contact_solution.at(dofs[local]);
    const ContactPointHistory& history = _contact_histories[constraint.contact][constraint.secondary];
    const Hex20AveragedConstraintValue value =
        hex20_averaged_constraint_value(constraint, state, committed_state, history);
    residual.assign(local_size, 0.0);
    if (jacobian != nullptr) jacobian->assign(local_size * local_size, 0.0);
    if (!(value.pressure > 0.0)) return;
    const std::array<double, 3> normal = {constraint.normal.x, constraint.normal.y, constraint.normal.z};
    for (std::size_t component = 0; component < 3; ++component) {
        for (std::size_t node = 0; node < node_count; ++node) {
            const std::size_t local = component * node_count + node;
            residual[local] = -constraint.area * constraint.gap_coefficients[node] *
                              (value.pressure * normal[component] + value.tangential_traction[component]);
        }
    }
    if (jacobian == nullptr) return;
    const NormalContactProperties& properties = _mechanical_properties[constraint.contact];
    std::array<double, 3> sliding_direction{};
    if (value.sliding)
        for (std::size_t component = 0; component < 3; ++component)
            sliding_direction[component] =
                value.trial_tangential_traction[component] / value.trial_tangential_magnitude;
    for (std::size_t row_component = 0; row_component < 3; ++row_component)
        for (std::size_t row_node = 0; row_node < node_count; ++row_node) {
            const std::size_t row = row_component * node_count + row_node;
            const double row_coefficient = constraint.gap_coefficients[row_node];
            for (std::size_t column_component = 0; column_component < 3; ++column_component)
                for (std::size_t column_node = 0; column_node < node_count; ++column_node) {
                    const std::size_t column = column_component * node_count + column_node;
                    const double column_coefficient = constraint.gap_coefficients[column_node],
                                 tangential_column_coefficient = -column_coefficient,
                                 pressure_derivative =
                                     -properties.penalty * normal[column_component] * column_coefficient;
                    double traction_derivative = normal[row_component] * pressure_derivative;
                    if (properties.friction_coefficient > 0.0) {
                        const double tangent_projector = (row_component == column_component ? 1.0 : 0.0) -
                                                         normal[row_component] * normal[column_component];
                        if (!value.sliding) {
                            traction_derivative +=
                                value.stick_stiffness * tangent_projector * tangential_column_coefficient;
                            if (properties.maximum_elastic_slip > 0.0)
                                traction_derivative +=
                                    value.tangential_traction[row_component] / value.pressure * pressure_derivative;
                        } else {
                            traction_derivative +=
                                properties.friction_coefficient *
                                (sliding_direction[row_component] * pressure_derivative +
                                    value.pressure * value.stick_stiffness / value.trial_tangential_magnitude *
                                        (tangent_projector -
                                            sliding_direction[row_component] * sliding_direction[column_component]) *
                                        tangential_column_coefficient);
                        }
                    }
                    (*jacobian)[row * local_size + column] = -constraint.area * row_coefficient * traction_derivative;
                }
        }
}

Quad8SurfaceContactLocalValues SpatialAssembly::hex20_contact_state(
    const Hex20ThermalCandidate& candidate, const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count()) throw std::invalid_argument("HEX20 contact state has the wrong size");
    const Quad8SurfaceContactLocalDofs dofs = hex20_contact_dofs(candidate);
    Quad8SurfaceContactLocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

Quad8SurfaceContactLocalValues SpatialAssembly::hex20_contact_state(
    const Hex20MechanicalCandidate& candidate, const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count()) throw std::invalid_argument("HEX20 contact state has the wrong size");
    const Quad8SurfaceContactLocalDofs dofs = hex20_contact_dofs(candidate);
    Quad8SurfaceContactLocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

void SpatialAssembly::build_contacts(const UnstructuredHex8Mesh& source_mesh) {
    _primary_boundaries.reserve(_definition.contacts.size());
    _secondary_boundaries.reserve(_definition.contacts.size());
    _primary_contact_faces.reserve(_definition.contacts.size());
    _secondary_contact_faces.reserve(_definition.contacts.size());
    _thermal_properties.reserve(_definition.contacts.size());
    _mechanical_properties.reserve(_definition.contacts.size());
    _thermal_point_counts.reserve(_definition.contacts.size());
    _contact_histories.resize(_definition.contacts.size());
    for (std::size_t contact_value = 0; contact_value < _definition.contacts.size(); ++contact_value) {
        ContactDefinition& definition = _definition.contacts[contact_value];
        if (definition.mechanical_discretization == MechanicalContactDiscretization::surface_to_surface)
            throw std::invalid_argument(
                "surface_to_surface mechanical contact currently requires HEX20 faces: " + definition.name);
        if (definition.friction_elastic_slip > 0.0)
            throw std::invalid_argument(
                "elastic_slip currently requires HEX20 surface_to_surface contact: " + definition.name);
        if (definition.mechanical_discretization == MechanicalContactDiscretization::automatic)
            definition.mechanical_discretization = MechanicalContactDiscretization::node_to_surface;
        if (definition.quad8_nodal_area_rule != Quad8NodalAreaRule::positive_lumped)
            throw std::invalid_argument(
                "quad8_nodal_area_rule = consistent_shape is supported only for HEX20 contact comparisons: " +
                definition.name);
        ResolvedBoundary primary = resolve_boundary(source_mesh, definition.primary),
                         secondary = resolve_boundary(source_mesh, definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("Three-dimensional self-contact is not supported: " + definition.name);
        const Hex8RegionMesh& primary_mesh = _meshes[primary.region];
        const Hex8RegionMesh& secondary_mesh = _meshes[secondary.region];
        for (const std::size_t secondary_local : secondary.boundary.nodes) {
            const std::size_t source_node = secondary_mesh.source_node_ids().at(secondary_local);
            if (std::any_of(
                    primary.boundary.nodes.begin(), primary.boundary.nodes.end(), [&](std::size_t primary_local) {
                        return primary_mesh.source_node_ids().at(primary_local) == source_node;
                    }))
                throw std::invalid_argument(
                    "Three-dimensional contact boundaries must not share source nodes: " + definition.name);
        }
        if (definition.mechanical_formulation == MechanicalContactFormulation::augmented_lagrangian)
            throw std::invalid_argument(
                "Three-dimensional contact currently supports the penalty formulation only: " + definition.name);
        if (definition.mechanical && definition.automatic_penalty) {
            const double primary_length = minimum_normal_length(primary_mesh, primary.boundary, definition.primary),
                         secondary_length =
                             minimum_normal_length(secondary_mesh, secondary.boundary, definition.secondary),
                         primary_modulus = _definition.regions[primary.region].material.reference_young_modulus,
                         secondary_modulus = _definition.regions[secondary.region].material.reference_young_modulus;
            definition.penalty =
                definition.penalty_factor / (primary_length / primary_modulus + secondary_length / secondary_modulus);
            if (!std::isfinite(definition.penalty) || !(definition.penalty > 0.0))
                throw std::overflow_error(
                    "Automatic three-dimensional contact penalty is not finite and positive: " + definition.name);
        }
        _thermal_properties.push_back({definition.thermal ? definition.gap_conductivity : 1.0,
            definition.thermal ? definition.minimum_gap : 1.0});
        _mechanical_properties.push_back({definition.mechanical ? definition.penalty : 1.0,
            definition.mechanical ? definition.friction_coefficient : 0.0, false});
        _contact_histories[contact_value].resize(secondary.boundary.nodes.size());
        std::vector<PrimaryContactFace> primary_faces;
        primary_faces.reserve(primary.boundary.faces.size());
        for (const Quad4FaceElement& primary_face : primary.boundary.faces) {
            const Quad4FaceCoordinates coordinates = face_coordinates(primary_mesh, primary_face);
            std::array<std::size_t, 4> nodes{};
            for (std::size_t node = 0; node < nodes.size(); ++node)
                nodes[node] = global_node(primary.region, primary_face.nodes[node]);
            primary_faces.push_back({nodes, coordinates, element_centroid(primary_mesh, primary_face.parent_element)});
        }
        std::vector<SecondaryContactFace> secondary_faces;
        secondary_faces.reserve(secondary.boundary.faces.size());
        std::size_t thermal_point_count = 0;
        for (const Quad4FaceElement& secondary_face : secondary.boundary.faces) {
            const Quad4FaceCoordinates secondary_coordinates = face_coordinates(secondary_mesh, secondary_face);
            const Quad4FaceGeometry secondary_geometry = make_quad4_face_geometry(secondary_coordinates);
            const CartesianPoint3 secondary_parent = element_centroid(secondary_mesh, secondary_face.parent_element);
            std::array<std::size_t, 4> secondary_nodes{};
            for (std::size_t node = 0; node < 4; ++node)
                secondary_nodes[node] = global_node(secondary.region, secondary_face.nodes[node]);
            const std::size_t secondary_face_index = secondary_faces.size();
            secondary_faces.push_back({secondary_nodes, secondary_coordinates, secondary_geometry, secondary_parent});
            if (definition.thermal) thermal_point_count += secondary_geometry.points.size();
            if (definition.mechanical)
                for (std::size_t secondary_local_node = 0; secondary_local_node < 4; ++secondary_local_node) {
                    const auto found = std::find(secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(),
                        secondary_face.nodes[secondary_local_node]);
                    if (found == secondary.boundary.nodes.end())
                        throw std::logic_error("Three-dimensional secondary contact node mapping failed");
                    const std::size_t secondary_node_index =
                        static_cast<std::size_t>(found - secondary.boundary.nodes.begin());
                    _mechanical_points.push_back(
                        {contact_value, secondary_node_index, secondary_face_index, secondary_local_node});
                }
        }
        _thermal_point_counts.push_back(thermal_point_count);
        _primary_contact_faces.push_back(std::move(primary_faces));
        _secondary_contact_faces.push_back(std::move(secondary_faces));
        _primary_boundaries.push_back(std::move(primary));
        _secondary_boundaries.push_back(std::move(secondary));
    }
    _thermal_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    _mechanical_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    _sparsity_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        _thermal_contact_offsets[contact + 1] = _thermal_contact_offsets[contact] + _thermal_point_counts[contact];
        _mechanical_contact_offsets[contact + 1] =
            _mechanical_contact_offsets[contact] + _contact_histories[contact].size();
        const std::size_t primary_count = _primary_contact_faces[contact].size(),
                          secondary_count = _secondary_contact_faces[contact].size();
        if (secondary_count > std::numeric_limits<std::size_t>::max() / primary_count)
            throw std::overflow_error("Three-dimensional contact sparsity candidate count exceeds size_t range");
        const std::size_t pairs = (_definition.contacts[contact].thermal || _definition.contacts[contact].mechanical)
                                      ? secondary_count * primary_count
                                      : 0;
        if (_sparsity_contact_offsets[contact] > std::numeric_limits<std::size_t>::max() - pairs)
            throw std::overflow_error("Three-dimensional contact sparsity candidate offset exceeds size_t range");
        _sparsity_contact_offsets[contact + 1] = _sparsity_contact_offsets[contact] + pairs;
    }
    _touched_thermal_points.resize(_thermal_contact_offsets.back());
    _thermal_minimum_distance.resize(_thermal_contact_offsets.back());
    _thermal_active_primary.resize(_thermal_contact_offsets.back());
    _thermal_cached_primary.assign(_thermal_contact_offsets.back(), std::numeric_limits<std::size_t>::max());
    _touched_mechanical_nodes.resize(_mechanical_contact_offsets.back());
    _mechanical_minimum_distance.resize(_mechanical_contact_offsets.back());
    _mechanical_selected_primary.resize(_mechanical_contact_offsets.back());
    _mechanical_cached_primary.assign(_mechanical_contact_offsets.back(), std::numeric_limits<std::size_t>::max());
    _mechanical_active_primary.resize(_mechanical_points.size());
    _contact_search_trees.resize(_definition.contacts.size());
}

void SpatialAssembly::build_hex20_contacts(const UnstructuredHex20Mesh& source_mesh) {
    _hex20_primary_boundaries.reserve(_definition.contacts.size());
    _hex20_secondary_boundaries.reserve(_definition.contacts.size());
    _hex20_primary_contact_faces.reserve(_definition.contacts.size());
    _hex20_secondary_contact_faces.reserve(_definition.contacts.size());
    _thermal_properties.reserve(_definition.contacts.size());
    _mechanical_properties.reserve(_definition.contacts.size());
    _thermal_point_counts.reserve(_definition.contacts.size());
    _contact_histories.resize(_definition.contacts.size());
    for (std::size_t contact_value = 0; contact_value < _definition.contacts.size(); ++contact_value) {
        ContactDefinition& definition = _definition.contacts[contact_value];
        ResolvedHex20Boundary primary = resolve_boundary(source_mesh, definition.primary),
                              secondary = resolve_boundary(source_mesh, definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("HEX20 self-contact is not supported: " + definition.name);
        const Hex20RegionMesh& primary_mesh = _hex20_meshes[primary.region];
        const Hex20RegionMesh& secondary_mesh = _hex20_meshes[secondary.region];
        for (const std::size_t secondary_local : secondary.boundary.displacement_nodes) {
            const std::size_t source_node = secondary_mesh.source_node_ids().at(secondary_local);
            if (std::any_of(primary.boundary.displacement_nodes.begin(), primary.boundary.displacement_nodes.end(),
                    [&](std::size_t primary_local) {
                        return primary_mesh.source_node_ids().at(primary_local) == source_node;
                    }))
                throw std::invalid_argument("HEX20 contact boundaries must not share source nodes: " + definition.name);
        }
        if (definition.mechanical_formulation == MechanicalContactFormulation::augmented_lagrangian)
            throw std::invalid_argument(
                "HEX20 contact currently supports the penalty formulation only: " + definition.name);
        if (definition.mechanical_discretization == MechanicalContactDiscretization::automatic)
            definition.mechanical_discretization = MechanicalContactDiscretization::surface_to_surface;
        const bool surface_to_surface =
            definition.mechanical_discretization == MechanicalContactDiscretization::surface_to_surface;
        // Abaqus-style averaged constraints represent small-strain,
        // small-sliding HEX20 contact.  Their relative vector displacement is
        // also the kinematic basis for the node-centered Coulomb law.
        const bool abaqus_averaged =
            surface_to_surface && _definition.regions[primary.region].strain_formulation == StrainFormulation::small &&
            _definition.regions[secondary.region].strain_formulation == StrainFormulation::small;
        if (definition.friction_elastic_slip > 0.0 && !surface_to_surface)
            throw std::invalid_argument("elastic_slip requires HEX20 surface_to_surface contact: " + definition.name);
        if (surface_to_surface && definition.quad8_nodal_area_rule != Quad8NodalAreaRule::positive_lumped)
            throw std::invalid_argument(
                "quad8_nodal_area_rule applies only to HEX20 node-to-surface contact comparisons: " + definition.name);
        if (definition.mechanical && definition.automatic_penalty) {
            const double primary_length = minimum_normal_length(primary_mesh, primary.boundary, definition.primary),
                         secondary_length =
                             minimum_normal_length(secondary_mesh, secondary.boundary, definition.secondary),
                         primary_modulus = _definition.regions[primary.region].material.reference_young_modulus,
                         secondary_modulus = _definition.regions[secondary.region].material.reference_young_modulus;
            definition.penalty =
                definition.penalty_factor / (primary_length / primary_modulus + secondary_length / secondary_modulus);
            if (!std::isfinite(definition.penalty) || !(definition.penalty > 0.0))
                throw std::overflow_error(
                    "Automatic HEX20 contact penalty is not finite and positive: " + definition.name);
        }
        _thermal_properties.push_back({definition.thermal ? definition.gap_conductivity : 1.0,
            definition.thermal ? definition.minimum_gap : 1.0});
        _mechanical_properties.push_back({definition.mechanical ? definition.penalty : 1.0,
            definition.mechanical ? definition.friction_coefficient : 0.0, false,
            definition.mechanical ? definition.friction_elastic_slip : 0.0});
        std::vector<Hex20PrimaryContactFace> primary_faces;
        primary_faces.reserve(primary.boundary.faces.size());
        for (const Quad8FaceElement& primary_face : primary.boundary.faces) {
            const Quad8FaceCoordinates coordinates = face_coordinates(primary_mesh, primary_face);
            std::array<std::size_t, 4> temperature_nodes{};
            std::array<std::size_t, 8> displacement_nodes{};
            for (std::size_t node = 0; node < 8; ++node) {
                displacement_nodes[node] = global_node(primary.region, primary_face.nodes[node]);
                if (node < 4)
                    temperature_nodes[node] = global_temperature_node(primary.region, primary_face.nodes[node]);
            }
            primary_faces.push_back({temperature_nodes, displacement_nodes, coordinates,
                hex20_element_centroid(primary_mesh, primary_face.parent_element)});
        }

        std::vector<Hex20SecondaryContactFace> secondary_faces;
        secondary_faces.reserve(secondary.boundary.faces.size());
        std::size_t thermal_point_count = 0, mechanical_point_count = 0;
        for (const Quad8FaceElement& secondary_face : secondary.boundary.faces) {
            const Quad8FaceCoordinates coordinates = face_coordinates(secondary_mesh, secondary_face);
            const Quad8FaceGeometry geometry = make_quad8_face_geometry(coordinates);
            std::array<std::size_t, 4> temperature_nodes{};
            std::array<std::size_t, 8> displacement_nodes{};
            for (std::size_t node = 0; node < 8; ++node) {
                displacement_nodes[node] = global_node(secondary.region, secondary_face.nodes[node]);
                if (node < 4)
                    temperature_nodes[node] = global_temperature_node(secondary.region, secondary_face.nodes[node]);
            }
            const std::size_t secondary_face_index = secondary_faces.size();
            secondary_faces.push_back({temperature_nodes, displacement_nodes, coordinates, geometry,
                hex20_element_centroid(secondary_mesh, secondary_face.parent_element), {}, {}, {}, {}, {}, {}});
            if (definition.thermal) thermal_point_count += geometry.thermal_points.size();
            if (definition.mechanical && !surface_to_surface)
                for (std::size_t secondary_local_node = 0; secondary_local_node < 8; ++secondary_local_node) {
                    const auto found = std::find(secondary.boundary.displacement_nodes.begin(),
                        secondary.boundary.displacement_nodes.end(), secondary_face.nodes[secondary_local_node]);
                    if (found == secondary.boundary.displacement_nodes.end())
                        throw std::logic_error("HEX20 secondary contact node mapping failed");
                    _hex20_mechanical_points.push_back(
                        {contact_value, static_cast<std::size_t>(found - secondary.boundary.displacement_nodes.begin()),
                            secondary_face_index, secondary_local_node, std::numeric_limits<std::size_t>::max()});
                }
        }
        if (definition.mechanical && surface_to_surface) {
            constexpr std::array<double, 3> points = {-0.7745966692414834, 0.0, 0.7745966692414834};
            constexpr std::array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

            struct Cell final {
                double xi_lower, xi_upper, eta_lower, eta_upper;
                std::size_t depth;
            };

            for (std::size_t secondary_face_index = 0; secondary_face_index < secondary_faces.size();
                ++secondary_face_index) {
                Hex20SecondaryContactFace& face = secondary_faces[secondary_face_index];
                const auto projection = [&](double xi, double eta, std::size_t primary_index) {
                    const Quad8FaceMechanicalQuadraturePoint point =
                        make_quad8_face_mechanical_point(face.coordinates, xi, eta, 1.0);
                    const Hex20PrimaryContactFace& primary_face = primary_faces.at(primary_index);
                    return compute_quad8_reference_projection(face.coordinates, primary_face.coordinates,
                        point.displacement_shape,
                        normal_orientation(
                            primary_face.coordinates, face.parent_centroid, primary_face.parent_centroid));
                };
                const auto owner = [&](double xi, double eta) {
                    const Quad8FaceMechanicalQuadraturePoint point =
                        make_quad8_face_mechanical_point(face.coordinates, xi, eta, 1.0);
                    double minimum_distance = std::numeric_limits<double>::infinity();
                    std::size_t selected = std::numeric_limits<std::size_t>::max();
                    for (std::size_t primary_index = 0; primary_index < primary_faces.size(); ++primary_index) {
                        const Hex20PrimaryContactFace& primary_face = primary_faces[primary_index];
                        const Quad8ReferenceProjectionValue reference = compute_quad8_reference_projection(
                            face.coordinates, primary_face.coordinates, point.displacement_shape,
                            normal_orientation(
                                primary_face.coordinates, face.parent_centroid, primary_face.parent_centroid));
                        if (!reference.projected) continue;
                        const double distance = std::abs(reference.gap);
                        if (distance < minimum_distance || (distance == minimum_distance && primary_index < selected)) {
                            minimum_distance = distance;
                            selected = primary_index;
                        }
                    }
                    return selected;
                };
                std::vector<Cell> pending = {{-1.0, 1.0, -1.0, 1.0, 0}}, leaves;
                while (!pending.empty()) {
                    const Cell cell = pending.back();
                    pending.pop_back();
                    const double xi_size = cell.xi_upper - cell.xi_lower, eta_size = cell.eta_upper - cell.eta_lower,
                                 xi_inset = 1.0e-9 * xi_size, eta_inset = 1.0e-9 * eta_size,
                                 xi_left = cell.xi_lower + xi_inset, xi_right = cell.xi_upper - xi_inset,
                                 eta_lower = cell.eta_lower + eta_inset, eta_upper = cell.eta_upper - eta_inset,
                                 xi_center = 0.5 * (cell.xi_lower + cell.xi_upper),
                                 eta_center = 0.5 * (cell.eta_lower + cell.eta_upper);
                    const std::array<std::size_t, 5> owners = {owner(xi_left, eta_lower), owner(xi_right, eta_lower),
                        owner(xi_right, eta_upper), owner(xi_left, eta_upper), owner(xi_center, eta_center)};
                    if (std::find(owners.begin(), owners.end(), std::numeric_limits<std::size_t>::max()) !=
                        owners.end())
                        throw std::invalid_argument("HEX20 small-sliding surface contact '" + definition.name +
                                                    "' has an unprojected reference averaging region");
                    const bool uniform = std::all_of(
                        owners.begin() + 1, owners.end(), [&](std::size_t value) { return value == owners.front(); });
                    if (uniform || cell.depth == 10) {
                        leaves.push_back(cell);
                        continue;
                    }
                    const bool corner_transition = owners[0] != owners[1] || owners[3] != owners[2] ||
                                                   owners[0] != owners[3] || owners[1] != owners[2],
                               split_xi = owners[0] != owners[1] || owners[3] != owners[2] || !corner_transition,
                               split_eta = owners[0] != owners[3] || owners[1] != owners[2] || !corner_transition;
                    if (split_xi && split_eta) {
                        pending.push_back({cell.xi_lower, xi_center, cell.eta_lower, eta_center, cell.depth + 1});
                        pending.push_back({xi_center, cell.xi_upper, cell.eta_lower, eta_center, cell.depth + 1});
                        pending.push_back({xi_center, cell.xi_upper, eta_center, cell.eta_upper, cell.depth + 1});
                        pending.push_back({cell.xi_lower, xi_center, eta_center, cell.eta_upper, cell.depth + 1});
                    } else if (split_xi) {
                        pending.push_back({cell.xi_lower, xi_center, cell.eta_lower, cell.eta_upper, cell.depth + 1});
                        pending.push_back({xi_center, cell.xi_upper, cell.eta_lower, cell.eta_upper, cell.depth + 1});
                    } else {
                        pending.push_back({cell.xi_lower, cell.xi_upper, cell.eta_lower, eta_center, cell.depth + 1});
                        pending.push_back({cell.xi_lower, cell.xi_upper, eta_center, cell.eta_upper, cell.depth + 1});
                    }
                }
                for (const Cell& cell : leaves) {
                    const double xi_half = 0.5 * (cell.xi_upper - cell.xi_lower),
                                 eta_half = 0.5 * (cell.eta_upper - cell.eta_lower),
                                 xi_center = 0.5 * (cell.xi_lower + cell.xi_upper),
                                 eta_center = 0.5 * (cell.eta_lower + cell.eta_upper);
                    for (std::size_t eta_point = 0; eta_point < points.size(); ++eta_point)
                        for (std::size_t xi_point = 0; xi_point < points.size(); ++xi_point) {
                            const double xi = xi_center + xi_half * points[xi_point],
                                         eta = eta_center + eta_half * points[eta_point],
                                         weight = xi_half * eta_half * weights[xi_point] * weights[eta_point];
                            const std::size_t primary_face = owner(xi, eta);
                            if (primary_face == std::numeric_limits<std::size_t>::max())
                                throw std::invalid_argument("HEX20 small-sliding surface contact '" + definition.name +
                                                            "' has an unprojected reference integration point");
                            face.contact_points.push_back(
                                make_quad8_face_mechanical_point(face.coordinates, xi, eta, weight));
                            face.contact_primary_faces.push_back(primary_face);
                            const Quad8ReferenceProjectionValue reference = projection(xi, eta, primary_face);
                            if (!reference.projected)
                                throw std::logic_error(
                                    "HEX20 small-sliding contact lost its reference projection while building");
                            face.contact_primary_shapes.push_back(reference.primary_shape);
                            face.contact_primary_derivatives_xi.push_back(reference.primary_derivative_xi);
                            face.contact_primary_derivatives_eta.push_back(reference.primary_derivative_eta);
                            face.contact_normal_orientations.push_back(
                                normal_orientation(primary_faces[primary_face].coordinates, face.parent_centroid,
                                    primary_faces[primary_face].parent_centroid));
                            if (!abaqus_averaged)
                                _hex20_mechanical_points.push_back({contact_value, mechanical_point_count,
                                    secondary_face_index, face.contact_points.size() - 1, primary_face});
                            ++mechanical_point_count;
                        }
                }
            }
        }
        if (definition.mechanical && abaqus_averaged) {
            struct ConstraintBuilder final {
                std::map<std::size_t, double> secondary, primary;
                std::map<std::size_t, CartesianPoint3> coordinates;
                std::map<std::size_t, double> secondary_output;
                CartesianPoint3 normal{};
                double area = 0.0;
            };

            std::vector<ConstraintBuilder> builders(secondary.boundary.displacement_nodes.size());
            const Matrix8& averaging = abaqus_quad8_averaging();
            for (std::size_t secondary_face_index = 0; secondary_face_index < secondary_faces.size();
                ++secondary_face_index) {
                const Hex20SecondaryContactFace& face = secondary_faces[secondary_face_index];
                const Quad8FaceElement& source_face = secondary.boundary.faces[secondary_face_index];
                Matrix8 gram{};
                double face_area = 0.0;
                for (const Quad8FaceMechanicalQuadraturePoint& point : face.geometry.mechanical_points) {
                    const double weighted_measure = point.quadrature_weight * reference_measure(point);
                    face_area += weighted_measure;
                    for (std::size_t row = 0; row < 8; ++row)
                        for (std::size_t column = 0; column < 8; ++column)
                            gram[row * 8 + column] +=
                                weighted_measure * point.displacement_shape[row] * point.displacement_shape[column];
                }
                if (!std::isfinite(face_area) || !(face_area > 0.0))
                    throw std::invalid_argument("HEX20 averaged contact has a nonpositive secondary face area");
                Matrix8 test_coefficients{};
                for (std::size_t local_constraint = 0; local_constraint < 8; ++local_constraint) {
                    const double fraction = local_constraint < 4 ? 1.0 / 24.0 : 5.0 / 24.0;
                    Vector8 weighted_averaging{};
                    for (std::size_t node = 0; node < 8; ++node)
                        weighted_averaging[node] = face_area * fraction * averaging[local_constraint * 8 + node];
                    const Vector8 coefficients = solve8(gram, weighted_averaging);
                    for (std::size_t node = 0; node < 8; ++node)
                        test_coefficients[local_constraint * 8 + node] = coefficients[node];

                    const auto found = std::find(secondary.boundary.displacement_nodes.begin(),
                        secondary.boundary.displacement_nodes.end(), source_face.nodes[local_constraint]);
                    if (found == secondary.boundary.displacement_nodes.end())
                        throw std::logic_error("HEX20 averaged contact secondary-node mapping failed");
                    const std::size_t output =
                        static_cast<std::size_t>(found - secondary.boundary.displacement_nodes.begin());
                    ConstraintBuilder& builder = builders[output];
                    const double local_area = face_area * fraction;
                    builder.area += local_area;
                    for (std::size_t node = 0; node < 8; ++node) {
                        const double value = weighted_averaging[node];
                        builder.secondary[face.displacement_nodes[node]] += value;
                        const auto output_node = std::find(secondary.boundary.displacement_nodes.begin(),
                            secondary.boundary.displacement_nodes.end(), source_face.nodes[node]);
                        if (output_node == secondary.boundary.displacement_nodes.end())
                            throw std::logic_error("HEX20 averaged contact output-node mapping failed");
                        builder.secondary_output[static_cast<std::size_t>(
                            output_node - secondary.boundary.displacement_nodes.begin())] += value;
                        builder.coordinates[face.displacement_nodes[node]] = face.coordinates[node];
                    }
                    const std::array<double, 2>& location = abaqus_quad8_constraint_locations()[local_constraint];
                    const Quad8FaceMechanicalQuadraturePoint constraint_point =
                        make_quad8_face_mechanical_point(face.coordinates, location[0], location[1], 1.0);
                    CartesianPoint3 face_normal = cross(constraint_point.tangent_xi, constraint_point.tangent_eta);
                    const std::size_t primary_index =
                        face.contact_primary_faces.at(face.contact_primary_faces.size() / 2);
                    const CartesianPoint3 direction =
                        subtract(primary_faces[primary_index].parent_centroid, face.parent_centroid);
                    if (dot(face_normal, direction) < 0.0) {
                        face_normal.x = -face_normal.x;
                        face_normal.y = -face_normal.y;
                        face_normal.z = -face_normal.z;
                    }
                    const double normal_measure = std::sqrt(dot(face_normal, face_normal));
                    builder.normal.x += local_area * face_normal.x / normal_measure;
                    builder.normal.y += local_area * face_normal.y / normal_measure;
                    builder.normal.z += local_area * face_normal.z / normal_measure;
                }

                for (std::size_t point_index = 0; point_index < face.contact_points.size(); ++point_index) {
                    const Quad8FaceMechanicalQuadraturePoint& point = face.contact_points[point_index];
                    const std::size_t primary_index = face.contact_primary_faces[point_index];
                    const Hex20PrimaryContactFace& primary_face = primary_faces[primary_index];
                    const double weighted_measure = point.quadrature_weight * reference_measure(point);
                    for (std::size_t local_constraint = 0; local_constraint < 8; ++local_constraint) {
                        double test = 0.0;
                        for (std::size_t node = 0; node < 8; ++node)
                            test += test_coefficients[local_constraint * 8 + node] * point.displacement_shape[node];
                        const auto found = std::find(secondary.boundary.displacement_nodes.begin(),
                            secondary.boundary.displacement_nodes.end(), source_face.nodes[local_constraint]);
                        const std::size_t output =
                            static_cast<std::size_t>(found - secondary.boundary.displacement_nodes.begin());
                        ConstraintBuilder& builder = builders[output];
                        for (std::size_t node = 0; node < 8; ++node) {
                            builder.primary[primary_face.displacement_nodes[node]] +=
                                weighted_measure * test * face.contact_primary_shapes[point_index][node];
                            builder.coordinates[primary_face.displacement_nodes[node]] = primary_face.coordinates[node];
                        }
                    }
                }
            }
            for (std::size_t output = 0; output < builders.size(); ++output) {
                const ConstraintBuilder& builder = builders[output];
                if (!std::isfinite(builder.area) || !(builder.area > 0.0))
                    throw std::invalid_argument("HEX20 averaged contact has a nonpositive constraint area");
                const double normal_measure = std::sqrt(dot(builder.normal, builder.normal));
                if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
                    throw std::invalid_argument("HEX20 averaged contact has an undefined constraint normal");
                Hex20AveragedConstraint constraint{};
                constraint.contact = contact_value;
                constraint.secondary = output;
                constraint.normal = {builder.normal.x / normal_measure, builder.normal.y / normal_measure,
                    builder.normal.z / normal_measure};
                constraint.area = builder.area;
                double secondary_sum = 0.0, primary_sum = 0.0;
                for (const auto& entry : builder.secondary) {
                    const double coefficient = -entry.second / builder.area;
                    if (std::abs(coefficient) <= 1.0e-14) continue;
                    constraint.nodes.push_back(entry.first);
                    constraint.gap_coefficients.push_back(coefficient);
                    secondary_sum -= coefficient;
                }
                for (const auto& entry : builder.primary) {
                    const double coefficient = entry.second / builder.area;
                    if (std::abs(coefficient) <= 1.0e-14) continue;
                    constraint.nodes.push_back(entry.first);
                    constraint.gap_coefficients.push_back(coefficient);
                    primary_sum += coefficient;
                }
                for (const auto& entry : builder.secondary_output) {
                    constraint.secondary_output_nodes.push_back(entry.first);
                    constraint.secondary_coefficients.push_back(entry.second / builder.area);
                }
                if (std::abs(secondary_sum - 1.0) > 1.0e-10 || std::abs(primary_sum - 1.0) > 1.0e-10)
                    throw std::invalid_argument("HEX20 averaged contact does not preserve rigid translation");
                CartesianPoint3 separation{};
                for (std::size_t node = 0; node < constraint.nodes.size(); ++node) {
                    const CartesianPoint3& coordinate = builder.coordinates.at(constraint.nodes[node]);
                    separation.x += constraint.gap_coefficients[node] * coordinate.x;
                    separation.y += constraint.gap_coefficients[node] * coordinate.y;
                    separation.z += constraint.gap_coefficients[node] * coordinate.z;
                }
                constraint.reference_gap = dot(separation, constraint.normal);
                _hex20_averaged_constraints.push_back(std::move(constraint));
            }
        }
        _contact_histories[contact_value].resize(
            !definition.mechanical
                ? 0
                : (surface_to_surface && !abaqus_averaged ? mechanical_point_count
                                                          : secondary.boundary.displacement_nodes.size()));
        _thermal_point_counts.push_back(thermal_point_count);
        _hex20_primary_contact_faces.push_back(std::move(primary_faces));
        _hex20_secondary_contact_faces.push_back(std::move(secondary_faces));
        _hex20_primary_boundaries.push_back(std::move(primary));
        _hex20_secondary_boundaries.push_back(std::move(secondary));
    }
    _thermal_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    _mechanical_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    _sparsity_contact_offsets.assign(_definition.contacts.size() + 1, 0);
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        _thermal_contact_offsets[contact + 1] = _thermal_contact_offsets[contact] + _thermal_point_counts[contact];
        _mechanical_contact_offsets[contact + 1] =
            _mechanical_contact_offsets[contact] + _contact_histories[contact].size();
        const std::size_t primary_count = _hex20_primary_contact_faces[contact].size(),
                          secondary_count = _hex20_secondary_contact_faces[contact].size();
        if (primary_count != 0 && secondary_count > std::numeric_limits<std::size_t>::max() / primary_count)
            throw std::overflow_error("HEX20 contact sparsity candidate count exceeds size_t range");
        const std::size_t pairs = (_definition.contacts[contact].thermal || _definition.contacts[contact].mechanical)
                                      ? secondary_count * primary_count
                                      : 0;
        if (_sparsity_contact_offsets[contact] > std::numeric_limits<std::size_t>::max() - pairs)
            throw std::overflow_error("HEX20 contact sparsity candidate offset exceeds size_t range");
        _sparsity_contact_offsets[contact + 1] = _sparsity_contact_offsets[contact] + pairs;
    }
    _touched_thermal_points.resize(_thermal_contact_offsets.back());
    _thermal_minimum_distance.resize(_thermal_contact_offsets.back());
    _thermal_active_primary.resize(_thermal_contact_offsets.back());
    _thermal_cached_primary.assign(_thermal_contact_offsets.back(), std::numeric_limits<std::size_t>::max());
    _touched_mechanical_nodes.resize(_mechanical_contact_offsets.back());
    _mechanical_minimum_distance.resize(_mechanical_contact_offsets.back());
    _mechanical_selected_primary.resize(_mechanical_contact_offsets.back());
    _mechanical_cached_primary.assign(_mechanical_contact_offsets.back(), std::numeric_limits<std::size_t>::max());
    _mechanical_active_primary.resize(_hex20_mechanical_points.size());
    _contact_search_trees.resize(_definition.contacts.size());
}

SpatialAssembly::ThermalCandidate SpatialAssembly::thermal_candidate(std::size_t point, std::size_t primary) const {
    const auto location =
        offset_location(_thermal_contact_offsets, point, "Three-dimensional thermal-contact point is out of range");
    const std::size_t contact = location.first, local_point = location.second, secondary_face = local_point / 4,
                      quadrature_point = local_point % 4;
    const SecondaryContactFace& secondary = _secondary_contact_faces.at(contact).at(secondary_face);
    const PrimaryContactFace& primary_face = _primary_contact_faces.at(contact).at(primary);
    std::array<std::size_t, 8> nodes{};
    std::copy(secondary.nodes.begin(), secondary.nodes.end(), nodes.begin());
    std::copy(primary_face.nodes.begin(), primary_face.nodes.end(), nodes.begin() + 4);
    const Quad4FaceQuadraturePoint& quadrature = secondary.geometry.points[quadrature_point];
    return {contact, nodes,
        {secondary.coordinates, primary_face.coordinates, quadrature.shape, quadrature.derivative_xi,
            quadrature.derivative_eta,
            normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid)},
        primary};
}

SpatialAssembly::MechanicalCandidate SpatialAssembly::mechanical_candidate(
    std::size_t point, std::size_t primary) const {
    const MechanicalPoint& metadata = _mechanical_points.at(point);
    const SecondaryContactFace& secondary = _secondary_contact_faces.at(metadata.contact).at(metadata.secondary_face);
    const PrimaryContactFace& primary_face = _primary_contact_faces.at(metadata.contact).at(primary);
    std::array<std::size_t, 8> nodes{};
    std::copy(secondary.nodes.begin(), secondary.nodes.end(), nodes.begin());
    std::copy(primary_face.nodes.begin(), primary_face.nodes.end(), nodes.begin() + 4);
    std::array<std::array<double, 4>, 4> shapes{}, derivatives_xi{}, derivatives_eta{};
    for (std::size_t q = 0; q < secondary.geometry.points.size(); ++q) {
        shapes[q] = secondary.geometry.points[q].shape;
        derivatives_xi[q] = secondary.geometry.points[q].derivative_xi;
        derivatives_eta[q] = secondary.geometry.points[q].derivative_eta;
    }
    return {metadata.contact, nodes,
        {secondary.coordinates, primary_face.coordinates, shapes, derivatives_xi, derivatives_eta,
            metadata.secondary_local_node,
            normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid)},
        metadata.secondary, primary};
}

SpatialAssembly::Hex20ThermalCandidate SpatialAssembly::hex20_thermal_candidate(
    std::size_t point, std::size_t primary) const {
    const auto location =
        offset_location(_thermal_contact_offsets, point, "HEX20 thermal-contact point is out of range");
    const std::size_t contact = location.first, local_point = location.second, secondary_face = local_point / 4,
                      quadrature_point = local_point % 4;
    const Hex20SecondaryContactFace& secondary = _hex20_secondary_contact_faces.at(contact).at(secondary_face);
    const Hex20PrimaryContactFace& primary_face = _hex20_primary_contact_faces.at(contact).at(primary);
    const Quad8FaceThermalQuadraturePoint& quadrature = secondary.geometry.thermal_points[quadrature_point];
    return {contact, secondary.temperature_nodes, primary_face.temperature_nodes, secondary.displacement_nodes,
        primary_face.displacement_nodes,
        {secondary.coordinates, primary_face.coordinates, quadrature.temperature_shape, quadrature.displacement_shape,
            quadrature.derivative_xi, quadrature.derivative_eta, quadrature.quadrature_weight,
            normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid)},
        primary};
}

SpatialAssembly::Hex20MechanicalCandidate SpatialAssembly::hex20_mechanical_candidate(
    std::size_t point, std::size_t primary) const {
    const Hex20MechanicalPoint& metadata = _hex20_mechanical_points.at(point);
    const Hex20SecondaryContactFace& secondary =
        _hex20_secondary_contact_faces.at(metadata.contact).at(metadata.secondary_face);
    const Hex20PrimaryContactFace& primary_face = _hex20_primary_contact_faces.at(metadata.contact).at(primary);
    const double orientation =
        normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid);
    Hex20MechanicalCandidate result{};
    result.contact = metadata.contact;
    result.secondary_temperature_nodes = secondary.temperature_nodes;
    result.primary_temperature_nodes = primary_face.temperature_nodes;
    result.secondary_displacement_nodes = secondary.displacement_nodes;
    result.primary_displacement_nodes = primary_face.displacement_nodes;
    result.surface_to_surface = _definition.contacts[metadata.contact].mechanical_discretization ==
                                MechanicalContactDiscretization::surface_to_surface;
    result.secondary = metadata.secondary;
    result.primary = primary;
    result.secondary_face = metadata.secondary_face;
    result.secondary_local_point = metadata.secondary_local_point;
    if (result.surface_to_surface) {
        const Quad8FaceMechanicalQuadraturePoint& quadrature =
            secondary.contact_points.at(metadata.secondary_local_point);
        result.surface_geometry = {secondary.coordinates, primary_face.coordinates, quadrature.displacement_shape,
            quadrature.derivative_xi, quadrature.derivative_eta,
            secondary.contact_primary_shapes.at(metadata.secondary_local_point),
            secondary.contact_primary_derivatives_xi.at(metadata.secondary_local_point),
            secondary.contact_primary_derivatives_eta.at(metadata.secondary_local_point), quadrature.quadrature_weight,
            secondary.contact_normal_orientations.at(metadata.secondary_local_point)};
        return result;
    }
    result.node_geometry = {secondary.coordinates, primary_face.coordinates, {}, {}, {}, {},
        metadata.secondary_local_point, orientation, _definition.contacts[metadata.contact].quad8_nodal_area_rule};
    for (std::size_t q = 0; q < quad8_surface_contact_quadrature_point_count; ++q) {
        const Quad8FaceMechanicalQuadraturePoint& quadrature = secondary.geometry.mechanical_points[q];
        result.node_geometry.secondary_shapes[q] = quadrature.displacement_shape;
        result.node_geometry.secondary_derivatives_xi[q] = quadrature.derivative_xi;
        result.node_geometry.secondary_derivatives_eta[q] = quadrature.derivative_eta;
        result.node_geometry.secondary_quadrature_weights[q] = quadrature.quadrature_weight;
    }
    return result;
}

SpatialAssembly::SparsityContact SpatialAssembly::sparsity_contact(std::size_t index) const {
    const auto location = offset_location(
        _sparsity_contact_offsets, index, "Three-dimensional sparsity contribution index is out of range");
    const std::size_t contact = location.first, primary_count = _primary_contact_faces.at(contact).size(),
                      secondary_face = location.second / primary_count, primary = location.second % primary_count;
    std::array<std::size_t, 8> nodes{};
    const std::array<std::size_t, 4>& secondary_nodes = _secondary_contact_faces.at(contact).at(secondary_face).nodes;
    const std::array<std::size_t, 4>& primary_nodes = _primary_contact_faces.at(contact).at(primary).nodes;
    std::copy(secondary_nodes.begin(), secondary_nodes.end(), nodes.begin());
    std::copy(primary_nodes.begin(), primary_nodes.end(), nodes.begin() + 4);
    return {nodes, _definition.contacts[contact].thermal, _definition.contacts[contact].mechanical};
}

SpatialAssembly::Hex20SparsityContact SpatialAssembly::hex20_sparsity_contact(std::size_t index) const {
    const auto location =
        offset_location(_sparsity_contact_offsets, index, "HEX20 sparsity contribution index is out of range");
    const std::size_t contact = location.first, primary_count = _hex20_primary_contact_faces.at(contact).size(),
                      secondary_face = location.second / primary_count, primary = location.second % primary_count;
    return {_hex20_secondary_contact_faces.at(contact).at(secondary_face).temperature_nodes,
        _hex20_primary_contact_faces.at(contact).at(primary).temperature_nodes,
        _hex20_secondary_contact_faces.at(contact).at(secondary_face).displacement_nodes,
        _hex20_primary_contact_faces.at(contact).at(primary).displacement_nodes, _definition.contacts[contact].thermal,
        _definition.contacts[contact].mechanical};
}

void SpatialAssembly::update_contact_search_trees(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional contact-search tree state size mismatch");
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        if (_uses_hex20) {
            const ResolvedHex20Boundary& primary = _hex20_primary_boundaries[contact];
            if (primary.boundary.faces.size() <= spatial_detail::contact_search_tree_minimum_items) continue;
            _contact_search_boxes.clear();
            _contact_search_boxes.reserve(primary.boundary.faces.size());
            for (std::size_t face_index = 0; face_index < primary.boundary.faces.size(); ++face_index) {
                const Quad8FaceElement& face = primary.boundary.faces[face_index];
                std::array<std::array<std::size_t, 3>, 8> displacement_dofs{};
                for (std::size_t node = 0; node < 8; ++node) {
                    const std::size_t global = global_node(primary.region, face.nodes[node]);
                    displacement_dofs[node] = {dof(Field::displacement_x, global), dof(Field::displacement_y, global),
                        dof(Field::displacement_z, global)};
                }
                spatial_detail::ContactSearchBox box =
                    quad8_search_box(face_coordinates(_hex20_meshes[primary.region], face), state, displacement_dofs);
                box.item = face_index;
                _contact_search_boxes.push_back(box);
            }
            if (_contact_search_trees[contact].can_refit(_contact_search_boxes.size()))
                _contact_search_trees[contact].refit(_contact_search_boxes);
            else
                _contact_search_trees[contact].build(_contact_search_boxes);
            continue;
        }
        const ResolvedBoundary& primary = _primary_boundaries[contact];
        if (primary.boundary.faces.size() <= spatial_detail::contact_search_tree_minimum_items) continue;
        const Hex8RegionMesh& mesh = _meshes[primary.region];
        _contact_search_boxes.clear();
        _contact_search_boxes.reserve(primary.boundary.faces.size());
        for (std::size_t face = 0; face < primary.boundary.faces.size(); ++face) {
            spatial_detail::ContactSearchBox box;
            box.minimum.fill(std::numeric_limits<double>::infinity());
            box.maximum.fill(-std::numeric_limits<double>::infinity());
            box.item = face;
            for (std::size_t local_node : primary.boundary.faces[face].nodes) {
                const std::size_t node = global_node(primary.region, local_node);
                const CartesianPoint3& reference = mesh.nodes().at(local_node);
                const std::array<double, 3> current = {reference.x + state[dof(Field::displacement_x, node)],
                    reference.y + state[dof(Field::displacement_y, node)],
                    reference.z + state[dof(Field::displacement_z, node)]};
                for (std::size_t component = 0; component < 3; ++component) {
                    box.minimum[component] = std::min(box.minimum[component], current[component]);
                    box.maximum[component] = std::max(box.maximum[component], current[component]);
                }
            }
            _contact_search_boxes.push_back(box);
        }
        if (_contact_search_trees[contact].can_refit(_contact_search_boxes.size()))
            _contact_search_trees[contact].refit(_contact_search_boxes);
        else
            _contact_search_trees[contact].build(_contact_search_boxes);
    }
}

std::size_t SpatialAssembly::mechanical_node_index(std::size_t contact, std::size_t node) const noexcept {
    return _mechanical_contact_offsets[contact] + node;
}

bool SpatialAssembly::mark_touched_thermal_points(std::size_t first, std::size_t last) const {
    std::fill(_touched_thermal_points.begin(), _touched_thermal_points.end(), 0U);
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t begin = std::max(first, ranges.thermal_begin), end = std::min(last, ranges.mechanical_begin);
    if (begin >= end) return false;
    for (std::size_t entry = begin; entry < end; ++entry) _touched_thermal_points[entry - ranges.thermal_begin] = 1U;
    return true;
}

bool SpatialAssembly::mark_touched_mechanical_nodes(std::size_t first, std::size_t last) const {
    std::fill(_touched_mechanical_nodes.begin(), _touched_mechanical_nodes.end(), 0U);
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t begin = std::max(first, ranges.mechanical_begin), end = std::min(last, ranges.boundary_begin);
    if (begin >= end) return false;
    if (_uses_hex20) {
        for (std::size_t entry = begin; entry < end; ++entry) {
            const std::size_t local = entry - ranges.mechanical_begin;
            if (local < _hex20_mechanical_points.size()) {
                const Hex20MechanicalPoint& point = _hex20_mechanical_points[local];
                _touched_mechanical_nodes[mechanical_node_index(point.contact, point.secondary)] = 1U;
            } else {
                const Hex20AveragedConstraint& constraint =
                    _hex20_averaged_constraints.at(local - _hex20_mechanical_points.size());
                _touched_mechanical_nodes[mechanical_node_index(constraint.contact, constraint.secondary)] = 1U;
            }
        }
    } else {
        for (std::size_t entry = begin; entry < end; ++entry) {
            const MechanicalPoint& point = _mechanical_points[entry - ranges.mechanical_begin];
            _touched_mechanical_nodes[mechanical_node_index(point.contact, point.secondary)] = 1U;
        }
    }
    return true;
}

void SpatialAssembly::update_thermal_candidates(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional thermal-contact shadow state size mismatch");
    if (!mark_touched_thermal_points(first, last)) return;
    if (_uses_hex20) {
        for (std::size_t point = 0; point < _thermal_active_primary.size(); ++point) {
            if (_touched_thermal_points[point] == 0U) continue;
            const auto location =
                offset_location(_thermal_contact_offsets, point, "HEX20 thermal-contact point is out of range");
            const std::size_t contact = location.first, primary_count = _hex20_primary_contact_faces[contact].size();
            _thermal_minimum_distance[point] = std::numeric_limits<double>::infinity();
            _thermal_active_primary[point] = std::numeric_limits<std::size_t>::max();
            const auto consider = [this, point, &state](std::size_t primary) {
                const Hex20ThermalCandidate candidate = hex20_thermal_candidate(point, primary);
                const ContactProjectionValue value =
                    compute_quad8_to_quad8_heat_projection(candidate.geometry, hex20_contact_state(candidate, state));
                if (!value.projected) return;
                const double distance = std::abs(value.gap);
                const std::size_t selected = _thermal_active_primary[point];
                if (distance < _thermal_minimum_distance[point] ||
                    (distance == _thermal_minimum_distance[point] &&
                        (selected == std::numeric_limits<std::size_t>::max() || primary < selected))) {
                    _thermal_minimum_distance[point] = distance;
                    _thermal_active_primary[point] = primary;
                }
            };
            const std::size_t cached_primary = _thermal_cached_primary[point];
            if (cached_primary < primary_count) consider(cached_primary);
            if (primary_count > spatial_detail::contact_search_tree_minimum_items) {
                const Hex20ThermalCandidate representative = hex20_thermal_candidate(point, 0);
                _contact_search_trees[contact].begin_query(
                    hex20_thermal_search_point(representative.geometry, hex20_contact_state(representative, state)),
                    _contact_search_query);
                std::size_t primary = 0;
                while (_contact_search_trees[contact].next_candidate(
                    _contact_search_query, _thermal_minimum_distance[point], primary)) {
                    if (primary != cached_primary) consider(primary);
                }
            } else {
                for (std::size_t primary = 0; primary < primary_count; ++primary)
                    if (primary != cached_primary) consider(primary);
            }
            if (_thermal_active_primary[point] != std::numeric_limits<std::size_t>::max())
                _thermal_cached_primary[point] = _thermal_active_primary[point];
        }
        return;
    }
    for (std::size_t point = 0; point < _thermal_active_primary.size(); ++point) {
        if (_touched_thermal_points[point] == 0U) continue;
        const auto location =
            offset_location(_thermal_contact_offsets, point, "Three-dimensional thermal-contact point is out of range");
        const std::size_t contact = location.first, primary_count = _primary_contact_faces[contact].size();
        _thermal_minimum_distance[point] = std::numeric_limits<double>::infinity();
        _thermal_active_primary[point] = std::numeric_limits<std::size_t>::max();
        const auto consider = [this, point, &state](std::size_t primary) {
            const ThermalCandidate candidate = thermal_candidate(point, primary);
            const ContactProjectionValue value =
                compute_quad4_to_quad4_heat_projection(candidate.geometry, contact_state(candidate.nodes, state));
            if (!value.projected) return;
            const double distance = std::abs(value.gap);
            const std::size_t selected = _thermal_active_primary[point];
            if (distance < _thermal_minimum_distance[point] ||
                (distance == _thermal_minimum_distance[point] &&
                    (selected == std::numeric_limits<std::size_t>::max() || primary < selected))) {
                _thermal_minimum_distance[point] = distance;
                _thermal_active_primary[point] = primary;
            }
        };
        const std::size_t cached_primary = _thermal_cached_primary[point];
        if (cached_primary < primary_count) consider(cached_primary);
        if (primary_count > spatial_detail::contact_search_tree_minimum_items) {
            const ThermalCandidate representative = thermal_candidate(point, 0);
            const Quad4SurfaceContactLocalValues representative_state = contact_state(representative.nodes, state);
            _contact_search_trees[contact].begin_query(
                thermal_search_point(representative.geometry, representative_state), _contact_search_query);
            std::size_t primary = 0;
            while (_contact_search_trees[contact].next_candidate(
                _contact_search_query, _thermal_minimum_distance[point], primary)) {
                if (primary != cached_primary) consider(primary);
            }
        } else
            for (std::size_t primary = 0; primary < primary_count; ++primary)
                if (primary != cached_primary) consider(primary);
        if (_thermal_active_primary[point] != std::numeric_limits<std::size_t>::max())
            _thermal_cached_primary[point] = _thermal_active_primary[point];
    }
}

void SpatialAssembly::update_mechanical_candidates(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional mechanical-contact shadow state size mismatch");
    if (!mark_touched_mechanical_nodes(first, last)) return;
    for (std::size_t node = 0; node < _mechanical_selected_primary.size(); ++node)
        if (_touched_mechanical_nodes[node] != 0U) {
            _mechanical_minimum_distance[node] = std::numeric_limits<double>::infinity();
            _mechanical_selected_primary[node] = std::numeric_limits<std::size_t>::max();
        }
    if (_uses_hex20) {
        for (std::size_t point = 0; point < _hex20_mechanical_points.size(); ++point) {
            const Hex20MechanicalPoint& metadata = _hex20_mechanical_points[point];
            const std::size_t node = mechanical_node_index(metadata.contact, metadata.secondary);
            _mechanical_active_primary[point] = std::numeric_limits<std::size_t>::max();
            if (_touched_mechanical_nodes[node] == 0U) continue;
            const bool surface_to_surface = _definition.contacts[metadata.contact].mechanical_discretization ==
                                            MechanicalContactDiscretization::surface_to_surface;
            if (surface_to_surface) {
                const Hex20MechanicalCandidate candidate =
                    hex20_mechanical_candidate(point, metadata.reference_primary);
                const ContactProjectionValue value = compute_quad8_to_quad8_contact_projection(
                    candidate.surface_geometry, hex20_contact_state(candidate, state));
                if (value.projected) {
                    _mechanical_minimum_distance[node] = std::abs(value.gap);
                    _mechanical_selected_primary[node] = metadata.reference_primary;
                    _mechanical_active_primary[point] = metadata.reference_primary;
                }
                continue;
            }
            const std::size_t primary_count = _hex20_primary_contact_faces[metadata.contact].size();
            const auto consider = [this, point, node, &state](std::size_t primary) {
                const Hex20MechanicalCandidate candidate = hex20_mechanical_candidate(point, primary);
                const Quad8SurfaceContactLocalValues candidate_state = hex20_contact_state(candidate, state);
                const ContactProjectionValue value =
                    candidate.surface_to_surface
                        ? compute_quad8_to_quad8_contact_projection(candidate.surface_geometry, candidate_state)
                        : compute_node_to_quad8_contact_projection(candidate.node_geometry, candidate_state);
                if (!value.projected) return;
                const double distance = std::abs(value.gap);
                if (distance < _mechanical_minimum_distance[node] ||
                    (distance == _mechanical_minimum_distance[node] && primary < _mechanical_selected_primary[node])) {
                    _mechanical_minimum_distance[node] = distance;
                    _mechanical_selected_primary[node] = primary;
                }
            };
            const std::size_t cached_primary = _mechanical_cached_primary[node];
            if (cached_primary < primary_count) consider(cached_primary);
            if (primary_count > spatial_detail::contact_search_tree_minimum_items) {
                const Hex20MechanicalCandidate representative = hex20_mechanical_candidate(point, 0);
                const Quad8SurfaceContactLocalValues representative_state = hex20_contact_state(representative, state);
                _contact_search_trees[metadata.contact].begin_query(
                    representative.surface_to_surface
                        ? hex20_mechanical_search_point(representative.surface_geometry, representative_state)
                        : hex20_mechanical_search_point(representative.node_geometry, representative_state),
                    _contact_search_query);
                std::size_t primary = 0;
                while (_contact_search_trees[metadata.contact].next_candidate(
                    _contact_search_query, _mechanical_minimum_distance[node], primary)) {
                    if (primary != cached_primary) consider(primary);
                }
            } else {
                for (std::size_t primary = 0; primary < primary_count; ++primary)
                    if (primary != cached_primary) consider(primary);
            }
        }
        for (std::size_t point = 0; point < _hex20_mechanical_points.size(); ++point) {
            const Hex20MechanicalPoint& metadata = _hex20_mechanical_points[point];
            const std::size_t node = mechanical_node_index(metadata.contact, metadata.secondary);
            const std::size_t primary = _mechanical_selected_primary[node];
            if (_touched_mechanical_nodes[node] == 0U || primary == std::numeric_limits<std::size_t>::max()) continue;
            const Hex20MechanicalCandidate candidate = hex20_mechanical_candidate(point, primary);
            const Quad8SurfaceContactLocalValues candidate_state = hex20_contact_state(candidate, state);
            const ContactProjectionValue value =
                candidate.surface_to_surface
                    ? compute_quad8_to_quad8_contact_projection(candidate.surface_geometry, candidate_state)
                    : compute_node_to_quad8_contact_projection(candidate.node_geometry, candidate_state);
            if (value.projected) _mechanical_active_primary[point] = primary;
        }
        for (const Hex20AveragedConstraint& constraint : _hex20_averaged_constraints) {
            const std::size_t node = mechanical_node_index(constraint.contact, constraint.secondary);
            if (_touched_mechanical_nodes[node] != 0U) _mechanical_selected_primary[node] = 0;
        }
        for (std::size_t node = 0; node < _mechanical_selected_primary.size(); ++node)
            if (_touched_mechanical_nodes[node] != 0U &&
                _mechanical_selected_primary[node] != std::numeric_limits<std::size_t>::max())
                _mechanical_cached_primary[node] = _mechanical_selected_primary[node];
        return;
    }
    for (std::size_t point = 0; point < _mechanical_points.size(); ++point) {
        const MechanicalPoint& metadata = _mechanical_points[point];
        const std::size_t node = mechanical_node_index(metadata.contact, metadata.secondary);
        _mechanical_active_primary[point] = std::numeric_limits<std::size_t>::max();
        if (_touched_mechanical_nodes[node] == 0U) continue;
        const std::size_t primary_count = _primary_contact_faces[metadata.contact].size();
        const auto consider = [this, point, node, &state](std::size_t primary) {
            const MechanicalCandidate candidate = mechanical_candidate(point, primary);
            const ContactProjectionValue value =
                compute_node_to_quad4_contact_projection(candidate.geometry, contact_state(candidate.nodes, state));
            if (!value.projected) return;
            const double distance = std::abs(value.gap);
            if (distance < _mechanical_minimum_distance[node] ||
                (distance == _mechanical_minimum_distance[node] && primary < _mechanical_selected_primary[node])) {
                _mechanical_minimum_distance[node] = distance;
                _mechanical_selected_primary[node] = primary;
            }
        };
        const std::size_t cached_primary = _mechanical_cached_primary[node];
        if (cached_primary < primary_count) consider(cached_primary);
        if (primary_count > spatial_detail::contact_search_tree_minimum_items) {
            const MechanicalCandidate representative = mechanical_candidate(point, 0);
            const Quad4SurfaceContactLocalValues representative_state = contact_state(representative.nodes, state);
            _contact_search_trees[metadata.contact].begin_query(
                mechanical_search_point(representative.geometry, representative_state), _contact_search_query);
            std::size_t primary = 0;
            while (_contact_search_trees[metadata.contact].next_candidate(
                _contact_search_query, _mechanical_minimum_distance[node], primary)) {
                if (primary != cached_primary) consider(primary);
            }
        } else
            for (std::size_t primary = 0; primary < primary_count; ++primary)
                if (primary != cached_primary) consider(primary);
    }
    for (std::size_t point = 0; point < _mechanical_points.size(); ++point) {
        const MechanicalPoint& metadata = _mechanical_points[point];
        const std::size_t node = mechanical_node_index(metadata.contact, metadata.secondary);
        const std::size_t primary = _mechanical_selected_primary[node];
        if (_touched_mechanical_nodes[node] == 0U || primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        const ContactProjectionValue value =
            compute_node_to_quad4_contact_projection(candidate.geometry, contact_state(candidate.nodes, state));
        if (value.projected) _mechanical_active_primary[point] = primary;
    }
    for (std::size_t node = 0; node < _mechanical_selected_primary.size(); ++node)
        if (_touched_mechanical_nodes[node] != 0U &&
            _mechanical_selected_primary[node] != std::numeric_limits<std::size_t>::max())
            _mechanical_cached_primary[node] = _mechanical_selected_primary[node];
}

std::vector<std::size_t> SpatialAssembly::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("Three-dimensional contribution range is out of bounds");
    std::vector<std::size_t> result;
    result.reserve(dof_count());
    std::vector<std::size_t> dofs;
    const auto append_face = [this, &result](const std::array<std::size_t, 4>& nodes) {
        for (const Field field :
            {Field::temperature, Field::displacement_x, Field::displacement_y, Field::displacement_z})
            for (const std::size_t node : nodes) result.push_back(dof(field, node));
    };
    const auto append_hex20_face = [this, &result](const std::array<std::size_t, 4>& temperature_nodes,
                                       const std::array<std::size_t, 8>& displacement_nodes) {
        for (const std::size_t node : temperature_nodes) result.push_back(dof(Field::temperature, node));
        for (const Field field : {Field::displacement_x, Field::displacement_y, Field::displacement_z})
            for (const std::size_t node : displacement_nodes) result.push_back(dof(field, node));
    };
    for (std::size_t entry = first; entry < last; ++entry) {
        contribution_dofs(entry, dofs);
        result.insert(result.end(), dofs.begin(), dofs.end());
    }
    std::vector<unsigned char> touched_contacts(_definition.contacts.size(), 0U);
    if (mark_touched_thermal_points(first, last)) {
        for (std::size_t point = 0; point < _touched_thermal_points.size(); ++point) {
            if (_touched_thermal_points[point] == 0U) continue;
            const auto location = offset_location(
                _thermal_contact_offsets, point, "Three-dimensional thermal-contact point is out of range");
            touched_contacts[location.first] = 1U;
            if (_uses_hex20) {
                const Hex20SecondaryContactFace& face =
                    _hex20_secondary_contact_faces[location.first][location.second / 4];
                append_hex20_face(face.temperature_nodes, face.displacement_nodes);
            } else
                append_face(_secondary_contact_faces[location.first][location.second / 4].nodes);
        }
    }
    if (mark_touched_mechanical_nodes(first, last)) {
        if (_uses_hex20) {
            for (const Hex20MechanicalPoint& point : _hex20_mechanical_points) {
                if (_touched_mechanical_nodes[mechanical_node_index(point.contact, point.secondary)] == 0U) continue;
                touched_contacts[point.contact] = 1U;
                const Hex20SecondaryContactFace& face =
                    _hex20_secondary_contact_faces[point.contact][point.secondary_face];
                append_hex20_face(face.temperature_nodes, face.displacement_nodes);
            }
        } else {
            for (const MechanicalPoint& point : _mechanical_points) {
                if (_touched_mechanical_nodes[mechanical_node_index(point.contact, point.secondary)] == 0U) continue;
                touched_contacts[point.contact] = 1U;
                append_face(_secondary_contact_faces[point.contact][point.secondary_face].nodes);
            }
        }
    }
    for (std::size_t contact = 0; contact < touched_contacts.size(); ++contact)
        if (touched_contacts[contact] != 0U) {
            if (_uses_hex20)
                for (const Hex20PrimaryContactFace& face : _hex20_primary_contact_faces[contact])
                    append_hex20_face(face.temperature_nodes, face.displacement_nodes);
            else
                for (const PrimaryContactFace& face : _primary_contact_faces[contact]) append_face(face.nodes);
        }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void SpatialAssembly::validate_local_state(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("Three-dimensional contribution range is out of bounds");
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional local validation state has the wrong size");
    update_contact_search_trees(state);
    update_thermal_candidates(first, last, state);
    update_mechanical_candidates(first, last, state);
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        if (_definition.contacts[contact].thermal) {
            std::size_t unprojected = 0;
            for (std::size_t point = _thermal_contact_offsets[contact]; point < _thermal_contact_offsets[contact + 1];
                ++point)
                if (_touched_thermal_points[point] != 0U &&
                    _thermal_active_primary[point] == std::numeric_limits<std::size_t>::max())
                    ++unprojected;
            if (unprojected != 0)
                throw std::domain_error("Three-dimensional thermal contact '" + _definition.contacts[contact].name +
                                        "' lost projection for " + std::to_string(unprojected) +
                                        " secondary integration points after searching every primary face");
        }
        if (_definition.contacts[contact].mechanical) {
            std::size_t unprojected = 0;
            for (std::size_t node = _mechanical_contact_offsets[contact];
                node < _mechanical_contact_offsets[contact + 1]; ++node)
                if (_touched_mechanical_nodes[node] != 0U &&
                    _mechanical_selected_primary[node] == std::numeric_limits<std::size_t>::max())
                    ++unprojected;
            if (unprojected != 0)
                throw std::domain_error("Three-dimensional mechanical contact '" + _definition.contacts[contact].name +
                                        "' lost projection for " + std::to_string(unprojected) +
                                        (_uses_hex20 && _definition.contacts[contact].mechanical_discretization ==
                                                            MechanicalContactDiscretization::surface_to_surface
                                                ? " fixed small-sliding integration points"
                                                : " secondary nodes after searching every primary face"));
        }
    }
}

void SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional committed contact state size mismatch");
    update_contact_search_trees(state);
    update_mechanical_candidates(0, contribution_count(), state);
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    std::vector<std::vector<bool>> updated(_definition.contacts.size());
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact)
        updated[contact].resize(_contact_histories[contact].size(), false);
    if (_uses_hex20) {
        for (std::size_t point = 0; point < _mechanical_active_primary.size(); ++point) {
            const std::size_t primary = _mechanical_active_primary[point];
            if (primary == std::numeric_limits<std::size_t>::max()) continue;
            const Hex20MechanicalCandidate candidate = hex20_mechanical_candidate(point, primary);
            const Quad8SurfaceContactLocalValues current = hex20_contact_state(candidate, state),
                                                 committed =
                                                     hex20_contact_state(candidate, _committed_contact_solution);
            const CartesianContactPointValue value =
                candidate.surface_to_surface
                    ? compute_quad8_to_quad8_contact_value(_mechanical_properties[candidate.contact],
                          candidate.surface_geometry, current, committed,
                          _contact_histories[candidate.contact][candidate.secondary])
                    : compute_node_to_quad8_contact_value(_mechanical_properties[candidate.contact],
                          candidate.node_geometry, current, committed,
                          _contact_histories[candidate.contact][candidate.secondary]);
            if (!value.projected) continue;
            ContactPointHistory trial = _contact_histories[candidate.contact][candidate.secondary];
            trial.sliding = value.sliding;
            trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
            if (updated[candidate.contact][candidate.secondary]) {
                const ContactPointHistory& prior = staged[candidate.contact][candidate.secondary];
                bool equal = prior.sliding == trial.sliding;
                for (std::size_t component = 0; component < 3; ++component) {
                    const double scale = std::max({1.0, std::abs(prior.cartesian_elastic_tangential_slip[component]),
                        std::abs(trial.cartesian_elastic_tangential_slip[component])});
                    equal = equal && std::abs(prior.cartesian_elastic_tangential_slip[component] -
                                              trial.cartesian_elastic_tangential_slip[component]) <=
                                         64.0 * std::numeric_limits<double>::epsilon() * scale;
                }
                if (!equal) throw std::logic_error("HEX20 secondary faces disagree on contact-node friction history");
                continue;
            }
            staged[candidate.contact][candidate.secondary] = trial;
            updated[candidate.contact][candidate.secondary] = true;
        }
        for (const Hex20AveragedConstraint& constraint : _hex20_averaged_constraints) {
            std::vector<std::size_t> dofs;
            hex20_averaged_constraint_dofs(constraint, dofs);
            std::vector<double> current(dofs.size()), committed(dofs.size());
            for (std::size_t local = 0; local < dofs.size(); ++local) {
                current[local] = state.at(dofs[local]);
                committed[local] = _committed_contact_solution.at(dofs[local]);
            }
            const Hex20AveragedConstraintValue value = hex20_averaged_constraint_value(
                constraint, current, committed, _contact_histories[constraint.contact][constraint.secondary]);
            ContactPointHistory trial = _contact_histories[constraint.contact][constraint.secondary];
            trial.sliding = value.sliding;
            trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
            if (updated[constraint.contact][constraint.secondary])
                throw std::logic_error("HEX20 averaged constraints share one friction-history slot");
            staged[constraint.contact][constraint.secondary] = trial;
            updated[constraint.contact][constraint.secondary] = true;
        }
        for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
            if (!_definition.contacts[contact].mechanical) continue;
            if (std::find(updated[contact].begin(), updated[contact].end(), false) != updated[contact].end())
                throw std::domain_error(
                    _definition.contacts[contact].mechanical_discretization ==
                            MechanicalContactDiscretization::surface_to_surface
                        ? "Cannot commit HEX20 friction history for an unprojected surface integration point"
                        : "Cannot commit HEX20 friction history for an unprojected node");
        }
        _contact_histories.swap(staged);
        _committed_contact_solution = state;
        return;
    }
    for (std::size_t point = 0; point < _mechanical_active_primary.size(); ++point) {
        const std::size_t primary = _mechanical_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        const CartesianContactPointValue value =
            compute_node_to_quad4_contact_value(_mechanical_properties[candidate.contact], candidate.geometry,
                contact_state(candidate.nodes, state), contact_state(candidate.nodes, _committed_contact_solution),
                _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        ContactPointHistory trial = _contact_histories[candidate.contact][candidate.secondary];
        trial.sliding = value.sliding;
        trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
        if (updated[candidate.contact][candidate.secondary]) {
            const ContactPointHistory& prior = staged[candidate.contact][candidate.secondary];
            bool equal = prior.sliding == trial.sliding;
            for (std::size_t component = 0; component < 3; ++component) {
                const double scale = std::max({1.0, std::abs(prior.cartesian_elastic_tangential_slip[component]),
                    std::abs(trial.cartesian_elastic_tangential_slip[component])});
                equal = equal && std::abs(prior.cartesian_elastic_tangential_slip[component] -
                                          trial.cartesian_elastic_tangential_slip[component]) <=
                                     64.0 * std::numeric_limits<double>::epsilon() * scale;
            }
            if (!equal)
                throw std::logic_error(
                    "Three-dimensional secondary half-faces disagree on contact-node friction history");
            continue;
        }
        staged[candidate.contact][candidate.secondary] = trial;
        updated[candidate.contact][candidate.secondary] = true;
    }
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        if (!_definition.contacts[contact].mechanical) continue;
        if (std::find(updated[contact].begin(), updated[contact].end(), false) != updated[contact].end())
            throw std::domain_error("Cannot commit three-dimensional friction history for an unprojected node");
    }
    _contact_histories.swap(staged);
    _committed_contact_solution = state;
}

void SpatialAssembly::restore_contact_state(
    const std::vector<double>& state, std::vector<std::vector<ContactPointHistory>> histories) {
    if (state.size() != dof_count() || histories.size() != _definition.contacts.size())
        throw std::invalid_argument("Three-dimensional restored contact state layout mismatch");
    for (std::size_t contact = 0; contact < histories.size(); ++contact) {
        if (histories[contact].size() != _contact_histories[contact].size())
            throw std::invalid_argument("Three-dimensional restored contact history layout mismatch");
        for (const ContactPointHistory& history : histories[contact]) {
            if (!std::isfinite(history.normal_multiplier) || history.normal_multiplier < 0.0)
                throw std::invalid_argument("Three-dimensional restored normal contact history is invalid");
            for (double component : history.cartesian_elastic_tangential_slip)
                if (!std::isfinite(component))
                    throw std::invalid_argument("Three-dimensional restored friction history is invalid");
        }
    }
    _committed_contact_solution = state;
    _contact_histories = std::move(histories);
}

std::vector<CartesianContactNodeSummary> SpatialAssembly::summarize_contact_nodes(
    std::size_t contact_value, const std::vector<double>& state) const {
    if (contact_value >= _definition.contacts.size())
        throw std::out_of_range("Three-dimensional contact index is out of range");
    update_contact_search_trees(state);
    update_mechanical_candidates(0, contribution_count(), state);
    if (_uses_hex20) {
        const ResolvedHex20Boundary& secondary = _hex20_secondary_boundaries.at(contact_value);
        const Hex20RegionMesh& mesh = _hex20_meshes.at(secondary.region);
        std::vector<CartesianContactNodeSummary> result;
        result.reserve(secondary.boundary.displacement_nodes.size());
        for (std::size_t node : secondary.boundary.displacement_nodes) {
            const CartesianPoint3& point = mesh.nodes().at(node);
            result.push_back({point.x, point.y, point.z, false, std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0, 0.0, 0.0, {}, {}, {}, {}, false});
        }
        const bool averaged = std::any_of(_hex20_averaged_constraints.begin(), _hex20_averaged_constraints.end(),
            [contact_value](const Hex20AveragedConstraint& value) { return value.contact == contact_value; });
        if (averaged) {
            std::vector<std::size_t> dofs;
            for (const Hex20AveragedConstraint& constraint : _hex20_averaged_constraints) {
                if (constraint.contact != contact_value) continue;
                hex20_averaged_constraint_dofs(constraint, dofs);
                std::vector<double> local_state(dofs.size());
                for (std::size_t local = 0; local < dofs.size(); ++local) local_state[local] = state.at(dofs[local]);
                std::vector<double> committed_state(dofs.size());
                for (std::size_t local = 0; local < dofs.size(); ++local)
                    committed_state[local] = _committed_contact_solution.at(dofs[local]);
                const Hex20AveragedConstraintValue value = hex20_averaged_constraint_value(constraint, local_state,
                    committed_state, _contact_histories[constraint.contact][constraint.secondary]);
                CartesianContactNodeSummary& predominant = result.at(constraint.secondary);
                predominant.projected = true;
                predominant.primary_face = 0;
                predominant.gap = value.gap;
                predominant.pressure = value.pressure;
                predominant.tributary_area = constraint.area;
                predominant.tangential_slip = value.tangential_slip;
                predominant.elastic_tangential_slip = value.elastic_tangential_slip;
                predominant.sliding = value.sliding;
                for (std::size_t entry = 0; entry < constraint.secondary_output_nodes.size(); ++entry) {
                    CartesianContactNodeSummary& output = result.at(constraint.secondary_output_nodes[entry]);
                    output.projected = true;
                    output.primary_face = 0;
                    output.contact_force += constraint.secondary_coefficients[entry] * value.force;
                    output.tangential_force += constraint.secondary_coefficients[entry] * value.tangential_force;
                    for (std::size_t component = 0; component < 3; ++component) {
                        const double normal = component == 0   ? constraint.normal.x
                                              : component == 1 ? constraint.normal.y
                                                               : constraint.normal.z;
                        output.normal_contact_force[component] +=
                            constraint.secondary_coefficients[entry] * value.force * normal;
                        output.tangential_contact_force[component] += constraint.secondary_coefficients[entry] *
                                                                      constraint.area *
                                                                      value.tangential_traction[component];
                    }
                }
            }
            for (CartesianContactNodeSummary& summary : result)
                if (summary.tributary_area > 0.0)
                    summary.tangential_traction = summary.tangential_force / summary.tributary_area;
            return result;
        }
        if (_definition.contacts[contact_value].mechanical_discretization ==
            MechanicalContactDiscretization::surface_to_surface) {
            std::vector<std::array<double, 3>> weighted_slip(result.size());
            std::vector<double> recovery_area(result.size()), recovery_force(result.size());
            for (std::size_t point = 0; point < _hex20_mechanical_points.size(); ++point) {
                const std::size_t primary = _mechanical_active_primary[point];
                if (primary == std::numeric_limits<std::size_t>::max()) continue;
                const Hex20MechanicalCandidate candidate = hex20_mechanical_candidate(point, primary);
                if (candidate.contact != contact_value) continue;
                const std::size_t history = mechanical_node_index(candidate.contact, candidate.secondary);
                if (candidate.primary != _mechanical_selected_primary[history]) continue;
                const Quad8SurfaceContactLocalValues current = hex20_contact_state(candidate, state),
                                                     committed =
                                                         hex20_contact_state(candidate, _committed_contact_solution);
                const CartesianContactPointValue value = compute_quad8_to_quad8_contact_value(
                    _mechanical_properties[candidate.contact], candidate.surface_geometry, current, committed,
                    _contact_histories[candidate.contact][candidate.secondary]);
                if (!value.projected) continue;
                const Quad8FaceElement& face = secondary.boundary.faces.at(candidate.secondary_face);
                for (std::size_t local_node = 0; local_node < 8; ++local_node) {
                    const auto found = std::find(secondary.boundary.displacement_nodes.begin(),
                        secondary.boundary.displacement_nodes.end(), face.nodes[local_node]);
                    if (found == secondary.boundary.displacement_nodes.end())
                        throw std::logic_error("HEX20 surface-contact output node mapping failed");
                    const std::size_t output_node =
                        static_cast<std::size_t>(found - secondary.boundary.displacement_nodes.begin());
                    CartesianContactNodeSummary& summary = result[output_node];
                    const double shape = candidate.surface_geometry.secondary_displacement_shape[local_node],
                                 nodal_area = shape * value.tributary_area;
                    summary.projected = true;
                    summary.primary_face = std::min(summary.primary_face, candidate.primary);
                    summary.gap = std::min(summary.gap, value.gap);
                    summary.tributary_area += nodal_area;
                    summary.contact_force += shape * value.contact_force;
                    for (std::size_t component = 0; component < 3; ++component)
                        summary.normal_contact_force[component] +=
                            shape * value.contact_force * value.normal[component];
                    summary.tangential_force += shape * value.tangential_force;
                    for (std::size_t component = 0; component < 3; ++component)
                        summary.tangential_contact_force[component] +=
                            shape * value.tributary_area * value.tangential_traction_vector[component];
                    for (std::size_t component = 0; component < 3; ++component)
                        weighted_slip[output_node][component] += nodal_area * value.elastic_tangential_slip[component];
                    summary.sliding = summary.sliding || value.sliding;
                }
            }
            NormalContactProperties recovery_properties = _mechanical_properties[contact_value];
            recovery_properties.friction_coefficient = 0.0;
            for (std::size_t secondary_face = 0; secondary_face < secondary.boundary.faces.size(); ++secondary_face) {
                const Hex20SecondaryContactFace& secondary_geometry =
                    _hex20_secondary_contact_faces[contact_value][secondary_face];
                const Quad8FaceElement& face = secondary.boundary.faces[secondary_face];
                for (const Quad8FaceMechanicalQuadraturePoint& quadrature :
                    secondary_geometry.geometry.mechanical_points) {
                    double minimum_distance = std::numeric_limits<double>::infinity();
                    std::size_t selected_primary = std::numeric_limits<std::size_t>::max();
                    Quad8ReferenceProjectionValue selected_projection{};
                    for (std::size_t primary_face = 0;
                        primary_face < _hex20_primary_contact_faces[contact_value].size(); ++primary_face) {
                        const Hex20PrimaryContactFace& primary_geometry =
                            _hex20_primary_contact_faces[contact_value][primary_face];
                        const Quad8ReferenceProjectionValue projection = compute_quad8_reference_projection(
                            secondary_geometry.coordinates, primary_geometry.coordinates, quadrature.displacement_shape,
                            normal_orientation(primary_geometry.coordinates, secondary_geometry.parent_centroid,
                                primary_geometry.parent_centroid));
                        if (!projection.projected || std::abs(projection.gap) >= minimum_distance) continue;
                        minimum_distance = std::abs(projection.gap);
                        selected_primary = primary_face;
                        selected_projection = projection;
                    }
                    if (selected_primary == std::numeric_limits<std::size_t>::max())
                        throw std::logic_error("HEX20 surface-contact pressure recovery has no primary face");
                    const Hex20PrimaryContactFace& primary_geometry =
                        _hex20_primary_contact_faces[contact_value][selected_primary];
                    Hex20MechanicalCandidate recovery{};
                    recovery.contact = contact_value;
                    recovery.secondary_temperature_nodes = secondary_geometry.temperature_nodes;
                    recovery.primary_temperature_nodes = primary_geometry.temperature_nodes;
                    recovery.secondary_displacement_nodes = secondary_geometry.displacement_nodes;
                    recovery.primary_displacement_nodes = primary_geometry.displacement_nodes;
                    recovery.surface_to_surface = true;
                    recovery.surface_geometry = {secondary_geometry.coordinates, primary_geometry.coordinates,
                        quadrature.displacement_shape, quadrature.derivative_xi, quadrature.derivative_eta,
                        selected_projection.primary_shape, selected_projection.primary_derivative_xi,
                        selected_projection.primary_derivative_eta, quadrature.quadrature_weight,
                        normal_orientation(primary_geometry.coordinates, secondary_geometry.parent_centroid,
                            primary_geometry.parent_centroid)};
                    const CartesianContactPointValue value = compute_quad8_to_quad8_contact_value(recovery_properties,
                        recovery.surface_geometry, hex20_contact_state(recovery, state),
                        hex20_contact_state(recovery, _committed_contact_solution), {});
                    for (std::size_t local_node = 0; local_node < face.nodes.size(); ++local_node) {
                        const auto found = std::find(secondary.boundary.displacement_nodes.begin(),
                            secondary.boundary.displacement_nodes.end(), face.nodes[local_node]);
                        if (found == secondary.boundary.displacement_nodes.end())
                            throw std::logic_error("HEX20 surface-contact pressure output node mapping failed");
                        const std::size_t output_node =
                            static_cast<std::size_t>(found - secondary.boundary.displacement_nodes.begin());
                        recovery_area[output_node] += quadrature.displacement_shape[local_node] * value.tributary_area;
                        recovery_force[output_node] += quadrature.displacement_shape[local_node] * value.contact_force;
                    }
                }
            }
            for (std::size_t node = 0; node < result.size(); ++node) {
                CartesianContactNodeSummary& summary = result[node];
                if (summary.tributary_area == 0.0) continue;
                if (recovery_area[node] != 0.0) summary.pressure = recovery_force[node] / recovery_area[node];
                summary.tangential_traction = summary.tangential_force / summary.tributary_area;
                for (std::size_t component = 0; component < 3; ++component)
                    summary.elastic_tangential_slip[component] =
                        weighted_slip[node][component] / summary.tributary_area;
            }
            return result;
        }
        for (std::size_t point = 0; point < _hex20_mechanical_points.size(); ++point) {
            const std::size_t primary = _mechanical_active_primary[point];
            if (primary == std::numeric_limits<std::size_t>::max()) continue;
            const Hex20MechanicalCandidate candidate = hex20_mechanical_candidate(point, primary);
            if (candidate.contact != contact_value) continue;
            const std::size_t node = mechanical_node_index(candidate.contact, candidate.secondary);
            if (candidate.primary != _mechanical_selected_primary[node]) continue;
            const CartesianContactPointValue value =
                compute_node_to_quad8_contact_value(_mechanical_properties[candidate.contact], candidate.node_geometry,
                    hex20_contact_state(candidate, state), hex20_contact_state(candidate, _committed_contact_solution),
                    _contact_histories[candidate.contact][candidate.secondary]);
            if (!value.projected) continue;
            CartesianContactNodeSummary& summary = result.at(candidate.secondary);
            if (summary.projected && summary.primary_face != candidate.primary)
                throw std::logic_error("HEX20 contact node has more than one active primary face");
            summary.projected = true;
            summary.primary_face = candidate.primary;
            summary.gap = std::min(summary.gap, value.gap);
            summary.tributary_area += value.tributary_area;
            summary.contact_force += value.contact_force;
            for (std::size_t component = 0; component < 3; ++component)
                summary.normal_contact_force[component] += value.contact_force * value.normal[component];
            summary.tangential_force += value.tangential_force;
            for (std::size_t component = 0; component < 3; ++component)
                summary.tangential_contact_force[component] +=
                    value.tributary_area * value.tangential_traction_vector[component];
            summary.elastic_tangential_slip = value.elastic_tangential_slip;
            summary.sliding = value.sliding;
        }
        for (CartesianContactNodeSummary& summary : result)
            if (summary.tributary_area != 0.0) {
                summary.pressure = summary.contact_force / summary.tributary_area;
                summary.tangential_traction = summary.tangential_force / summary.tributary_area;
            }
        return result;
    }
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const Hex8RegionMesh& mesh = _meshes.at(secondary.region);
    std::vector<CartesianContactNodeSummary> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) {
        const CartesianPoint3& point = mesh.nodes().at(node);
        result.push_back({point.x, point.y, point.z, false, std::numeric_limits<std::size_t>::max(),
            std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0, 0.0, 0.0, {}, {}, {}, {}, false});
    }
    for (std::size_t point = 0; point < _mechanical_active_primary.size(); ++point) {
        const std::size_t primary = _mechanical_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        if (candidate.contact != contact_value) continue;
        const std::size_t node = mechanical_node_index(candidate.contact, candidate.secondary);
        if (candidate.primary != _mechanical_selected_primary[node]) continue;
        const CartesianContactPointValue value =
            compute_node_to_quad4_contact_value(_mechanical_properties[candidate.contact], candidate.geometry,
                contact_state(candidate.nodes, state), contact_state(candidate.nodes, _committed_contact_solution),
                _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        CartesianContactNodeSummary& summary = result.at(candidate.secondary);
        if (summary.projected && summary.primary_face != candidate.primary)
            throw std::logic_error("Three-dimensional contact node has more than one active primary face");
        summary.projected = true;
        summary.primary_face = candidate.primary;
        summary.gap = std::min(summary.gap, value.gap);
        summary.tributary_area += value.tributary_area;
        summary.contact_force += value.contact_force;
        for (std::size_t component = 0; component < 3; ++component)
            summary.normal_contact_force[component] += value.contact_force * value.normal[component];
        summary.tangential_force += value.tangential_force;
        for (std::size_t component = 0; component < 3; ++component)
            summary.tangential_contact_force[component] +=
                value.tributary_area * value.tangential_traction_vector[component];
        summary.elastic_tangential_slip = value.elastic_tangential_slip;
        summary.sliding = value.sliding;
    }
    for (CartesianContactNodeSummary& summary : result)
        if (summary.tributary_area > 0.0) {
            summary.pressure = summary.contact_force / summary.tributary_area;
            summary.tangential_traction = summary.tangential_force / summary.tributary_area;
        }
    return result;
}

std::vector<std::size_t> SpatialAssembly::contact_secondary_source_nodes(std::size_t contact_value) const {
    if (_uses_hex20) {
        const ResolvedHex20Boundary& secondary = _hex20_secondary_boundaries.at(contact_value);
        const Hex20RegionMesh& mesh = _hex20_meshes.at(secondary.region);
        std::vector<std::size_t> result;
        result.reserve(secondary.boundary.displacement_nodes.size());
        for (std::size_t node : secondary.boundary.displacement_nodes)
            result.push_back(mesh.source_node_ids().at(node));
        return result;
    }
    const ResolvedBoundary& secondary = _secondary_boundaries.at(contact_value);
    const Hex8RegionMesh& mesh = _meshes.at(secondary.region);
    std::vector<std::size_t> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) result.push_back(mesh.source_node_ids().at(node));
    return result;
}

InterfaceSummary SpatialAssembly::summarize_interface(
    std::size_t contact_value, const std::vector<double>& state) const {
    if (contact_value >= _definition.contacts.size())
        throw std::out_of_range("Three-dimensional contact index is out of range");
    validate_state(state);
    InterfaceSummary result;
    for (std::size_t point = _thermal_contact_offsets[contact_value];
        point < _thermal_contact_offsets[contact_value + 1]; ++point) {
        const std::size_t primary = _thermal_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        if (_uses_hex20) {
            const Hex20ThermalCandidate candidate = hex20_thermal_candidate(point, primary);
            const CartesianHeatQuadratureValue value = compute_quad8_to_quad8_gap_heat_value(
                _thermal_properties[contact_value], candidate.geometry, hex20_contact_state(candidate, state));
            if (!value.projected) throw std::logic_error("Active HEX20 thermal-contact candidate is not projected");
            result.minimum_gap = std::min(result.minimum_gap, value.gap);
            result.total_heat_rate += value.weighted_measure * value.heat_flux;
            continue;
        }
        const ThermalCandidate candidate = thermal_candidate(point, primary);
        const CartesianHeatQuadratureValue value = compute_quad4_to_quad4_gap_heat_value(
            _thermal_properties[contact_value], candidate.geometry, contact_state(candidate.nodes, state));
        if (!value.projected)
            throw std::logic_error("Active three-dimensional thermal-contact candidate is not projected");
        result.minimum_gap = std::min(result.minimum_gap, value.gap);
        result.total_heat_rate += value.weighted_measure * value.heat_flux;
    }
    if (!_definition.contacts[contact_value].thermal) result.minimum_gap = 0.0;
    if (_definition.contacts[contact_value].mechanical) {
        for (const CartesianContactNodeSummary& node : summarize_contact_nodes(contact_value, state)) {
            if (!node.projected) {
                ++result.unprojected_contact_nodes;
                continue;
            }
            ++result.projected_contact_nodes;
            result.minimum_contact_gap = std::min(result.minimum_contact_gap, node.gap);
            result.maximum_contact_pressure = std::max(result.maximum_contact_pressure, node.pressure);
            result.total_contact_force += node.contact_force;
            result.total_tangential_force += node.tangential_force;
            if (node.pressure > 0.0) ++result.active_contact_nodes;
        }
    } else {
        result.minimum_contact_gap = 0.0;
    }
    return result;
}

void SpatialAssembly::refresh_controls() {
    for (std::size_t region = 0; region < region_count(); ++region)
        _kernel_data[region].volumetric_heat_source = region_heat_source(region);
    refresh_dirichlet_values();
    for (std::size_t kernel = 0; kernel < _boundary_data.size(); ++kernel) {
        const BoundaryConditionDefinition& boundary =
            _definition.boundary_conditions[_boundary_definition_indices[kernel]];
        if (boundary.type != BoundaryConditionType::convection) {
            _boundary_data[kernel].load = spatial_detail::controlled_value(
                _definition, _time, _load_factor, boundary.value, boundary.scale_with_load, boundary.function);
            continue;
        }
        const spatial_detail::ConvectionValues values = convection_values(boundary);
        _boundary_data[kernel].load = values.coefficient;
        _boundary_data[kernel].ambient_temperature = values.ambient;
    }
}
} // namespace fuelsim::cartesian
