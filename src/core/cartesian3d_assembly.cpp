#include "cartesian3d_assembly.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
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

struct SurfacePullbackCovectors final {
    CartesianPoint3 first, second;
};

SurfacePullbackCovectors surface_pullback_covectors(const Quad4FaceQuadraturePoint& current,
    const Quad4FaceQuadraturePoint& reference, double normal_orientation, double tangent_orientation,
    const CartesianPoint3& current_normal, const CartesianPoint3& current_first) {
    const CartesianPoint3 current_second = cross(current_normal, current_first),
                          reference_area = cross(reference.tangent_xi, reference.tangent_eta);
    const double reference_area_measure = std::sqrt(dot(reference_area, reference_area)),
                 reference_xi_measure = std::sqrt(dot(reference.tangent_xi, reference.tangent_xi));
    if (!std::isfinite(reference_area_measure) || !(reference_area_measure > 0.0) ||
        !std::isfinite(reference_xi_measure) || !(reference_xi_measure > 0.0))
        throw std::domain_error("Three-dimensional averaged friction has a degenerate reference tangent metric");
    const CartesianPoint3 reference_normal = {normal_orientation * reference_area.x / reference_area_measure,
                              normal_orientation * reference_area.y / reference_area_measure,
                              normal_orientation * reference_area.z / reference_area_measure},
                          reference_first = {tangent_orientation * reference.tangent_xi.x / reference_xi_measure,
                              tangent_orientation * reference.tangent_xi.y / reference_xi_measure,
                              tangent_orientation * reference.tangent_xi.z / reference_xi_measure};
    const CartesianPoint3 reference_second = cross(reference_normal, reference_first);
    const double current_11 = dot(current.tangent_xi, current_first),
                 current_12 = dot(current.tangent_eta, current_first),
                 current_21 = dot(current.tangent_xi, current_second),
                 current_22 = dot(current.tangent_eta, current_second),
                 reference_11 = dot(reference.tangent_xi, reference_first),
                 reference_12 = dot(reference.tangent_eta, reference_first),
                 reference_21 = dot(reference.tangent_xi, reference_second),
                 reference_22 = dot(reference.tangent_eta, reference_second),
                 current_determinant = current_11 * current_22 - current_12 * current_21;
    if (!std::isfinite(current_determinant) || !(std::abs(current_determinant) > 0.0))
        throw std::domain_error("Three-dimensional averaged friction has a singular current tangent metric");
    const double first_first = (reference_11 * current_22 - reference_12 * current_21) / current_determinant,
                 first_second = (-reference_11 * current_12 + reference_12 * current_11) / current_determinant,
                 second_first = (reference_21 * current_22 - reference_22 * current_21) / current_determinant,
                 second_second = (-reference_21 * current_12 + reference_22 * current_11) / current_determinant;
    return {{first_first * current_first.x + first_second * current_second.x,
                first_first * current_first.y + first_second * current_second.y,
                first_first * current_first.z + first_second * current_second.z},
        {second_first * current_first.x + second_second * current_second.x,
            second_first * current_first.y + second_second * current_second.y,
            second_first * current_first.z + second_second * current_second.z}};
}

using Matrix4 = std::array<double, 16>;
using Vector4 = std::array<double, 4>;
using Matrix8 = std::array<double, 64>;

struct AbaqusQuad8TransferSample final {
    double first;
    double second;
    double weight;
};

const std::array<std::array<double, 2>, 4>& abaqus_quad4_constraint_locations() {
    // B3.8 identifies these effective centers for Abaqus/Standard C3D8
    // small-sliding surface-to-surface constraints in face-node order.
    static const std::array<std::array<double, 2>, 4> value = {
        {{{-0.5, -0.5}}, {{0.5, -0.5}}, {{0.5, 0.5}}, {{-0.5, 0.5}}}};
    return value;
}

const Matrix4& abaqus_quad4_averaging() {
    static const Matrix4 value = {{9.0 / 16.0, 3.0 / 16.0, 1.0 / 16.0, 3.0 / 16.0, 3.0 / 16.0, 9.0 / 16.0, 3.0 / 16.0,
        1.0 / 16.0, 1.0 / 16.0, 3.0 / 16.0, 9.0 / 16.0, 3.0 / 16.0, 3.0 / 16.0, 1.0 / 16.0, 3.0 / 16.0, 9.0 / 16.0}};
    return value;
}

Vector4 solve4(Matrix4 matrix, Vector4 right_hand_side) {
    for (std::size_t column = 0; column < 4; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < 4; ++row)
            if (std::abs(matrix[row * 4 + column]) > std::abs(matrix[pivot * 4 + column])) pivot = row;
        if (!(std::abs(matrix[pivot * 4 + column]) > 1.0e-14))
            throw std::invalid_argument("HEX8 averaged contact has a singular secondary face Gram matrix");
        if (pivot != column) {
            for (std::size_t entry = column; entry < 4; ++entry)
                std::swap(matrix[column * 4 + entry], matrix[pivot * 4 + entry]);
            std::swap(right_hand_side[column], right_hand_side[pivot]);
        }
        const double diagonal = matrix[column * 4 + column];
        for (std::size_t row = column + 1; row < 4; ++row) {
            const double factor = matrix[row * 4 + column] / diagonal;
            for (std::size_t entry = column; entry < 4; ++entry)
                matrix[row * 4 + entry] -= factor * matrix[column * 4 + entry];
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    Vector4 result{};
    for (std::size_t reverse = 0; reverse < 4; ++reverse) {
        const std::size_t row = 3 - reverse;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < 4; ++column) value -= matrix[row * 4 + column] * result[column];
        result[row] = value / matrix[row * 4 + row];
    }
    return result;
}

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

const Matrix8& abaqus_quad8_pressure_recovery() {
    // H20.42 identifies this as the exact linear stage in Abaqus/Standard: project the eight nonnegative constraint
    // pressures on each quadratic face onto the Q4-compatible bilinear nodal subspace.  The matrix is the equal-norm
    // least-squares projector.  The probe outputs retain additional state-dependent nonlinear behavior.
    static const Matrix8 value = {{17.0 / 24.0, -3.0 / 24.0, 1.0 / 24.0, -3.0 / 24.0, 7.0 / 24.0, -1.0 / 24.0,
        -1.0 / 24.0, 7.0 / 24.0, -3.0 / 24.0, 17.0 / 24.0, -3.0 / 24.0, 1.0 / 24.0, 7.0 / 24.0, 7.0 / 24.0, -1.0 / 24.0,
        -1.0 / 24.0, 1.0 / 24.0, -3.0 / 24.0, 17.0 / 24.0, -3.0 / 24.0, -1.0 / 24.0, 7.0 / 24.0, 7.0 / 24.0,
        -1.0 / 24.0, -3.0 / 24.0, 1.0 / 24.0, -3.0 / 24.0, 17.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, 7.0 / 24.0,
        7.0 / 24.0, 7.0 / 24.0, 7.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, 7.0 / 24.0, 3.0 / 24.0, -1.0 / 24.0, 3.0 / 24.0,
        -1.0 / 24.0, 7.0 / 24.0, 7.0 / 24.0, -1.0 / 24.0, 3.0 / 24.0, 7.0 / 24.0, 3.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0,
        -1.0 / 24.0, 7.0 / 24.0, 7.0 / 24.0, -1.0 / 24.0, 3.0 / 24.0, 7.0 / 24.0, 3.0 / 24.0, 7.0 / 24.0, -1.0 / 24.0,
        -1.0 / 24.0, 7.0 / 24.0, 3.0 / 24.0, -1.0 / 24.0, 3.0 / 24.0, 7.0 / 24.0}};
    return value;
}

const std::array<AbaqusQuad8TransferSample, 12>& abaqus_quad8_corner_transfer_rule() {
    // H20.28 identifies this mesh-independent parent-face sampling rule from
    // Abaqus/Standard R2018x primary meshes with 8, 16, and 32 subdivisions.
    // The first two coordinates use [0, 1] on the complete secondary face.
    static const std::array<AbaqusQuad8TransferSample, 12> value = {
        {{0.028276366456420152, 0.12450141132135768, 0.086805555555562158},
            {0.03199788301796376, 0.21899929433932117, 0.13765653331648353},
            {0.033360052620092687, 0.033360052620092881, 0.12916616234614192},
            {0.047248941508981852, 0.38603095325083253, 0.056145392295365942},
            {0.10552883626879452, 0.1055288362687966, 0.058333837653884016},
            {0.11941772515768481, 0.16952460230472302, 0.086215718815705361},
            {0.12450141132135777, 0.028276366456420156, 0.086805555555567654},
            {0.16952460230472319, 0.11941772515768576, 0.086215718815703918},
            {0.17633545031537054, 0.26711181677179063, 0.03942680001684136},
            {0.21899929433932119, 0.031997883017963316, 0.13765653331647812},
            {0.26711181677178519, 0.1763354503153691, 0.039426800016858818},
            {0.38603095325083064, 0.047248941508978792, 0.056145392295423285}}};
    return value;
}

const std::array<AbaqusQuad8TransferSample, 40>& abaqus_quad8_edge_transfer_rule() {
    static const std::array<AbaqusQuad8TransferSample, 40> value = {
        {{0.028276366456418445, 0.12450141132135609, 0.0014138183228190539},
            {0.031997883017964308, 0.21899929433932022, 0.0035378857178991023},
            {0.033360052620099695, 0.033360052620091216, 0.0019890098870065402},
            {0.047248941508978626, 0.38603095325083892, 0.0079520568710641521},
            {0.10552883626880082, 0.10552883626879425, 0.0042609901129871917},
            {0.10566243270259684, 0.5, 0.014174682452723866},
            {0.11941772515768409, 0.16952460230472471, 0.0097562764622762116},
            {0.12450141132135764, 0.028276366456420281, 0.0062250705660643595},
            {0.16952460230472319, 0.11941772515768448, 0.013849944679972161},
            {0.17633545031536924, 0.26711181677179047, 0.024587114282105537},
            {0.21899929433932128, 0.031997883017963434, 0.024213929253937019},
            {0.25, 0.35566243270259301, 0.052900635094618217}, {0.25, 0.64433756729741143, 0.0095993649054004881},
            {0.26711181677179036, 0.17633545031537012, 0.037244404079395804},
            {0.35566243270259157, 0.75000000000000011, 0.006250000000009805},
            {0.35566243270259346, 0.25, 0.081249999999999295},
            {0.38603095325083253, 0.047248941508981623, 0.064969499764469363},
            {0.39433756729740682, 0.5, 0.03582531754730961}, {0.5, 0.10566243270259353, 0.11997595264190444},
            {0.5, 0.39433756729740604, 0.055024047358088829}, {0.5, 0.60566243270259357, 0.023325317547314144},
            {0.5, 0.89433756729742642, 0.0016746824527028559}, {0.60566243270259401, 0.5, 0.03582531754730503},
            {0.61396904674916786, 0.047248941508981374, 0.064969499764467697},
            {0.64433756729740654, 0.25, 0.081249999999991149}, {0.64433756729740688, 0.75, 0.0062499999999983133},
            {0.73288818322821081, 0.17633545031536949, 0.037244404079388477},
            {0.74999999999999989, 0.35566243270259401, 0.052900635094603674},
            {0.75, 0.64433756729740743, 0.0095993649053841227},
            {0.78100070566067914, 0.03199788301796349, 0.024213929253929546},
            {0.82366454968462943, 0.26711181677178969, 0.024587114282093269},
            {0.83047539769527723, 0.11941772515768355, 0.013849944679969167},
            {0.8754985886786425, 0.028276366456420572, 0.0062250705660644168},
            {0.88058227484231988, 0.16952460230472366, 0.0097562764622718401},
            {0.89433756729740699, 0.5, 0.0141746824526763},
            {0.89447116373120827, 0.10552883626879543, 0.0042609901129861093},
            {0.95275105849101915, 0.38603095325082687, 0.0079520568710483106},
            {0.96663994737990877, 0.03336005262008996, 0.0019890098870046931},
            {0.96800211698203353, 0.21899929433932186, 0.0035378857178827283},
            {0.97172363354357416, 0.12450141132135775, 0.0014138183228165758}}};
    return value;
}

std::vector<AbaqusQuad8TransferSample> abaqus_quad8_primary_transfer_rule(std::size_t local_constraint) {
    std::vector<AbaqusQuad8TransferSample> result;
    if (local_constraint < 4) {
        result.reserve(abaqus_quad8_corner_transfer_rule().size());
        for (const AbaqusQuad8TransferSample& sample : abaqus_quad8_corner_transfer_rule()) result.push_back(sample);
    } else if (local_constraint < 8) {
        result.reserve(abaqus_quad8_edge_transfer_rule().size());
        for (const AbaqusQuad8TransferSample& sample : abaqus_quad8_edge_transfer_rule()) result.push_back(sample);
    } else {
        throw std::out_of_range("Abaqus HEX20 contact constraint index exceeds the Quad8 face");
    }
    for (AbaqusQuad8TransferSample& sample : result) {
        const double first = sample.first, second = sample.second;
        if (local_constraint == 1)
            sample.first = 1.0 - first;
        else if (local_constraint == 2) {
            sample.first = 1.0 - first;
            sample.second = 1.0 - second;
        } else if (local_constraint == 3)
            sample.second = 1.0 - second;
        else if (local_constraint == 5) {
            sample.first = 1.0 - second;
            sample.second = first;
        } else if (local_constraint == 6)
            sample.second = 1.0 - second;
        else if (local_constraint == 7) {
            sample.first = second;
            sample.second = first;
        }
    }
    const double sum = std::accumulate(result.begin(), result.end(), 0.0,
        [](double value, const AbaqusQuad8TransferSample& sample) { return value + sample.weight; });
    if (!std::isfinite(sum) || !(sum > 0.0))
        throw std::logic_error("Abaqus HEX20 primary transfer rule has an invalid weight sum");
    for (AbaqusQuad8TransferSample& sample : result) sample.weight /= sum;
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

std::array<double, 3> mechanical_search_point(
    const Quad4ToQuad4MechanicalGeometry& geometry, const Quad4SurfaceContactLocalValues& state) {
    std::array<double, 3> result{};
    for (std::size_t node = 0; node < 4; ++node) {
        result[0] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].x + state[8 + node]);
        result[1] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].y + state[16 + node]);
        result[2] += geometry.secondary_shape[node] * (geometry.secondary_coordinates[node].z + state[24 + node]);
    }
    return result;
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

double primary_material_orientation(
    const Quad4FaceCoordinates& primary_coordinates, const CartesianPoint3& primary_parent_centroid) {
    return normal_orientation(primary_coordinates, face_centroid(primary_coordinates), primary_parent_centroid);
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

double primary_material_orientation(
    const Quad8FaceCoordinates& primary_coordinates, const CartesianPoint3& primary_parent_centroid) {
    return normal_orientation(primary_coordinates, quad8_face_centroid(primary_coordinates), primary_parent_centroid);
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
    if (boundary.type == BoundaryConditionType::heat_flux)
        return {Quad4FaceBoundaryKind::surface_heat_flux, CartesianTractionComponent::x, boundary.value, 0.0,
            use_displaced_geometry};
    return {Quad4FaceBoundaryKind::convection, CartesianTractionComponent::x, boundary.heat_transfer_coefficient,
        boundary.ambient_temperature, use_displaced_geometry};
}

SpatialContributionType boundary_contribution_type(BoundaryConditionType type) {
    if (type == BoundaryConditionType::pressure) return SpatialContributionType::pressure;
    if (type == BoundaryConditionType::traction) return SpatialContributionType::traction;
    if (type == BoundaryConditionType::heat_flux) return SpatialContributionType::heat_flux;
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
        if (boundary.type == BoundaryConditionType::pressure || boundary.type == BoundaryConditionType::traction ||
            boundary.type == BoundaryConditionType::heat_flux || boundary.type == BoundaryConditionType::convection)
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
            region_heat_source(region_index), 0.0, region(region_index).strain_formulation,
            region(region_index).hex8_element_formulation, region(region_index).initial_temperature});
    }
    _committed_contact_solution = initial_state();
    validate_local_state(0, contribution_count(), _committed_contact_solution);
    refresh_controls();
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh)
    : SpatialLayout(definition, spatial_detail::resolve_block_ids(definition, source_mesh, true, true),
          spatial_detail::DofLayout::cartesian_3d),
      _uses_hex20(true) {
    for (const RegionDefinition& region : _definition.regions)
        if (region.hex8_element_formulation != Hex8ElementFormulation::c3d8t)
            throw std::invalid_argument("C3D8RT formulation requires an eight-node HEX8 mesh: " + region.name);
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
        if (boundary.type == BoundaryConditionType::heat_flux)
            throw std::invalid_argument("HEX20 surface heat flux is not implemented");
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
        if (boundary.type == BoundaryConditionType::convection && boundary.configuration_explicit)
            throw std::invalid_argument(
                "Convection configuration is currently supported only for C3D8T and C3D8RT: " + boundary.name);
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
            region_heat_source(region), 0.0, this->region(region).strain_formulation,
            this->region(region).hex8_element_formulation, this->region(region).initial_temperature});
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

void SpatialAssembly::set_heat_source_interval(double begin_time, double end_time) {
    for (std::size_t region = 0; region < region_count(); ++region)
        _kernel_data[region].volumetric_heat_source = region_heat_source_average(region, begin_time, end_time);
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
        const std::size_t point_count = _uses_hex20 ? _hex20_mechanical_points.size() : _mechanical_points.size();
        if (point >= point_count) {
            averaged_constraint_dofs(_abaqus_averaged_constraints.at(point - point_count), dofs);
            return;
        }
        if (_uses_hex20) {
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
        const auto location = element_location(index);
        const RegionDefinition& element_region = region(location.first);
        const bool finite_reduced = element_region.hex8_element_formulation == Hex8ElementFormulation::c3d8rt &&
                                    element_region.strain_formulation == StrainFormulation::finite;
        const std::size_t thermal_column_end = finite_reduced ? hex8_local_dof_count : 8;
        set_pattern_block(pattern, hex8_local_dof_count, 0, 8, 0, thermal_column_end);
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
        const std::size_t point = index - ranges.mechanical_begin,
                          point_count = _uses_hex20 ? _hex20_mechanical_points.size() : _mechanical_points.size();
        if (point >= point_count) {
            const AbaqusAveragedConstraint& constraint = _abaqus_averaged_constraints.at(point - point_count);
            const std::size_t local_size =
                3 * (constraint.active_nodes.empty() ? constraint.nodes.size() : constraint.active_nodes.size());
            pattern.assign(local_size * local_size, 1U);
            return;
        }
        if (_uses_hex20) {
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
        if (data.use_displaced_geometry) set_pattern_block(pattern, quad4_face_local_dof_count, 0, 4, 4, 16);
    } else if (data.kind == Quad4FaceBoundaryKind::surface_heat_flux && data.use_displaced_geometry) {
        set_pattern_block(pattern, quad4_face_local_dof_count, 0, 4, 4, 16);
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
    if (index >= _sparsity_contact_offsets.back()) {
        averaged_sparsity_contribution_dofs(index - _sparsity_contact_offsets.back(), dofs);
        return;
    }
    if (_uses_hex20) {
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
    if (index >= _sparsity_contact_offsets.back()) {
        std::vector<std::size_t> dofs;
        averaged_sparsity_contribution_dofs(index - _sparsity_contact_offsets.back(), dofs);
        const std::size_t local_size = dofs.size();
        pattern.assign(local_size * local_size, 1U);
        return;
    }
    if (_uses_hex20) {
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
        const std::size_t point = index - ranges.mechanical_begin,
                          point_count = _uses_hex20 ? _hex20_mechanical_points.size() : _mechanical_points.size();
        if (point >= point_count) {
            compute_averaged_constraint(
                _abaqus_averaged_constraints.at(point - point_count), state, residual, jacobian);
            return;
        }
        if (_uses_hex20) {
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
        const MechanicalCandidate entry = mechanical_candidate(point, _mechanical_active_primary.at(point));
        Quad4SurfaceContactLocalValues current{};
        std::copy(state.begin(), state.end(), current.begin());
        const Quad4SurfaceContactLocalValues committed = contribution_state(index, _committed_contact_solution);
        NormalContactProperties point_properties = _mechanical_properties[entry.contact];
        if (entry.surface_to_surface &&
            std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
                [&entry](const AbaqusAveragedConstraint& constraint) {
                    return constraint.contact == entry.contact && constraint.friction_only;
                })) {
            point_properties.friction_coefficient = 0.0;
            point_properties.maximum_elastic_slip = 0.0;
        }
        Quad4SurfaceContactLocalJacobian local_jacobian{};
        const Quad4SurfaceContactLocalResidual result =
            entry.surface_to_surface
                ? compute_quad4_to_quad4_contact(point_properties, entry.surface_geometry, current, committed,
                      _contact_histories[entry.contact][entry.secondary],
                      jacobian == nullptr ? nullptr : &local_jacobian)
                : compute_node_to_quad4_contact(_mechanical_properties[entry.contact], entry.node_geometry, current,
                      committed, _contact_histories[entry.contact][entry.secondary],
                      jacobian == nullptr ? nullptr : &local_jacobian);
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

double SpatialAssembly::mechanical_hourglass_energy(
    std::size_t region, std::size_t element, const Hex8LocalValues& state) const {
    return compute_hex8_mechanical_hourglass_energy(
        _kernel_data.at(region), region_element_geometry(region, element), state);
}

SpatialAssembly::ContributionRanges SpatialAssembly::contribution_ranges() const noexcept {
    const std::size_t thermal_begin = volume_contribution_count(),
                      mechanical_begin = thermal_begin + _thermal_contact_offsets.back(),
                      mechanical_count = (_uses_hex20 ? _hex20_mechanical_points.size() : _mechanical_points.size()) +
                                         _abaqus_averaged_constraints.size(),
                      boundary_begin = mechanical_begin + mechanical_count;
    const std::size_t boundary_count =
        _uses_hex20 ? _hex20_boundary_contributions.size() : _boundary_contributions.size();
    return {thermal_begin, mechanical_begin, boundary_begin, boundary_begin + boundary_count};
}

std::size_t SpatialAssembly::sparsity_contribution_count() const noexcept {
    // Every face boundary block is a subset of its adjacent volume block. Contact quadrature points and face nodes
    // that share one secondary-face/primary-face pair also have the same 32-DOF graph, so one representative preserves
    // the complete graph without repeating it for each runtime contribution.
    return volume_contribution_count() + _sparsity_contact_offsets.back() + averaged_sparsity_contribution_count();
}

bool SpatialAssembly::jacobian_sparsity_is_state_dependent() const noexcept {
    return !_uses_hex20 && std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
                               [](const AbaqusAveragedConstraint& constraint) { return constraint.finite_sliding; });
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

void SpatialAssembly::averaged_constraint_dofs(
    const AbaqusAveragedConstraint& constraint, std::vector<std::size_t>& dofs) const {
    const std::vector<std::size_t>& nodes =
        constraint.active_nodes.empty() ? constraint.nodes : constraint.active_nodes;
    dofs.clear();
    dofs.reserve(3 * nodes.size());
    for (const Field field : {Field::displacement_x, Field::displacement_y, Field::displacement_z})
        for (const std::size_t node : nodes) dofs.push_back(dof(field, node));
}

std::size_t SpatialAssembly::averaged_sparsity_contribution_count() const noexcept {
    if (!_uses_hex20) {
        std::size_t result = 0;
        for (auto constraint = _abaqus_averaged_constraints.begin(); constraint != _abaqus_averaged_constraints.end();
            ++constraint) {
            if (!constraint->finite_sliding) {
                ++result;
                continue;
            }
            const bool first_for_contact = std::none_of(_abaqus_averaged_constraints.begin(), constraint,
                [constraint](const AbaqusAveragedConstraint& previous) {
                    return previous.finite_sliding && previous.contact == constraint->contact;
                });
            if (first_for_contact) ++result;
        }
        return result;
    }
    std::size_t result = _abaqus_averaged_constraints.size();
    for (auto constraint = _abaqus_averaged_constraints.begin(); constraint != _abaqus_averaged_constraints.end();
        ++constraint)
        if (constraint->finite_sliding && std::none_of(_abaqus_averaged_constraints.begin(), constraint,
                                              [constraint](const AbaqusAveragedConstraint& previous) {
                                                  return previous.finite_sliding &&
                                                         previous.contact == constraint->contact;
                                              }))
            ++result;
    return result;
}

void SpatialAssembly::averaged_sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    if (!_uses_hex20) {
        for (auto constraint = _abaqus_averaged_constraints.begin(); constraint != _abaqus_averaged_constraints.end();
            ++constraint) {
            const bool first_for_contact =
                !constraint->finite_sliding ||
                std::none_of(_abaqus_averaged_constraints.begin(), constraint,
                    [constraint](const AbaqusAveragedConstraint& previous) {
                        return previous.finite_sliding && previous.contact == constraint->contact;
                    });
            if (!first_for_contact) continue;
            if (index != 0) {
                --index;
                continue;
            }
            std::vector<std::size_t> nodes;
            if (constraint->finite_sliding) {
                for (const AbaqusAveragedConstraint& candidate : _abaqus_averaged_constraints)
                    if (candidate.finite_sliding && candidate.contact == constraint->contact)
                        nodes.insert(nodes.end(), candidate.nodes.begin(), candidate.nodes.end());
            } else {
                nodes = constraint->nodes;
            }
            std::sort(nodes.begin(), nodes.end());
            nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
            dofs.clear();
            dofs.reserve(3 * nodes.size());
            for (const Field field : {Field::displacement_x, Field::displacement_y, Field::displacement_z})
                for (const std::size_t node : nodes) dofs.push_back(dof(field, node));
            return;
        }
        throw std::out_of_range("Three-dimensional averaged-contact sparsity index is out of range");
    }
    for (auto constraint = _abaqus_averaged_constraints.begin(); constraint != _abaqus_averaged_constraints.end();
        ++constraint) {
        std::vector<std::size_t> nodes;
        if (constraint->finite_sliding) {
            for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint->finite_sliding_samples) {
                const Hex20SecondaryContactFace& secondary =
                    _hex20_secondary_contact_faces.at(constraint->contact).at(sample.secondary_face);
                nodes.insert(nodes.end(), secondary.displacement_nodes.begin(), secondary.displacement_nodes.end());
            }
        } else {
            nodes = constraint->nodes;
        }
        if (index == 0) {
            std::sort(nodes.begin(), nodes.end());
            nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
            dofs.clear();
            dofs.reserve(3 * nodes.size());
            for (const Field field : {Field::displacement_x, Field::displacement_y, Field::displacement_z})
                for (const std::size_t node : nodes) dofs.push_back(dof(field, node));
            return;
        }
        --index;
        const bool first_for_contact =
            constraint->finite_sliding && std::none_of(_abaqus_averaged_constraints.begin(), constraint,
                                              [constraint](const AbaqusAveragedConstraint& previous) {
                                                  return previous.finite_sliding &&
                                                         previous.contact == constraint->contact;
                                              });
        if (!first_for_contact) continue;
        if (index != 0) {
            --index;
            continue;
        }
        nodes.clear();
        for (const Hex20PrimaryContactFace& primary : _hex20_primary_contact_faces.at(constraint->contact))
            nodes.insert(nodes.end(), primary.displacement_nodes.begin(), primary.displacement_nodes.end());
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        dofs.clear();
        dofs.reserve(3 * nodes.size());
        for (const Field field : {Field::displacement_x, Field::displacement_y, Field::displacement_z})
            for (const std::size_t node : nodes) dofs.push_back(dof(field, node));
        return;
    }
    throw std::out_of_range("HEX20 averaged-contact sparsity index is out of range");
}

SpatialAssembly::AbaqusAveragedConstraintValue SpatialAssembly::averaged_constraint_value(
    const AbaqusAveragedConstraint& constraint, const std::vector<double>& state,
    const std::vector<double>& committed_state, const ContactPointHistory& history) const {
    if (constraint.finite_region_normal) return finite_region_normal_value(constraint, constraint.active_nodes, state);
    const std::size_t node_count =
        constraint.active_nodes.empty() ? constraint.nodes.size() : constraint.active_nodes.size();
    if (state.size() != 3 * node_count || committed_state.size() != state.size())
        throw std::invalid_argument("Abaqus-style averaged-contact state has the wrong size");
    const auto stored_node = [&constraint](std::size_t local) {
        return constraint.active_node_indices.empty() ? local : constraint.active_node_indices.at(local);
    };
    const bool local_normal_gradients = !constraint.normal_gap_coefficients.empty();
    if (local_normal_gradients && constraint.normal_gap_coefficients.size() != constraint.nodes.size())
        throw std::logic_error("Abaqus-style averaged contact has inconsistent local-normal gap coefficients");
    if (!constraint.projected) return {};
    std::array<double, 3> relative{}, committed_relative{};
    for (std::size_t component = 0; component < 3; ++component) {
        for (std::size_t node = 0; node < node_count; ++node) {
            const double coefficient = constraint.gap_coefficients[stored_node(node)];
            relative[component] += coefficient * state[component * node_count + node];
            committed_relative[component] += coefficient * committed_state[component * node_count + node];
        }
    }
    const std::array<double, 3> normal = {constraint.normal.x, constraint.normal.y, constraint.normal.z};
    AbaqusAveragedConstraintValue result{};
    result.tangent_first = {constraint.tangent_first.x, constraint.tangent_first.y, constraint.tangent_first.z};
    result.gap = constraint.reference_gap;
    if (local_normal_gradients)
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < node_count; ++node)
                result.gap += constraint.normal_gap_coefficients[stored_node(node)][component] *
                              state[component * node_count + node];
    else
        for (std::size_t component = 0; component < 3; ++component)
            result.gap += normal[component] * relative[component];
    double normal_relative = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        normal_relative += normal[component] * relative[component];
    for (std::size_t component = 0; component < 3; ++component)
        result.tangential_slip[component] = -relative[component] + normal_relative * normal[component];
    if (constraint.finite_sliding) {
        const std::array<double, 3> current_first = result.tangent_first,
                                    current_second = {normal[1] * current_first[2] - normal[2] * current_first[1],
                                        normal[2] * current_first[0] - normal[0] * current_first[2],
                                        normal[0] * current_first[1] - normal[1] * current_first[0]};
        double total_first = 0.0, total_second = 0.0;
        for (std::size_t component = 0; component < 3; ++component) {
            total_first += history.cartesian_total_tangential_slip[component] * current_first[component];
            total_second += history.cartesian_total_tangential_slip[component] * current_second[component];
        }
        for (std::size_t component = 0; component < 3; ++component)
            result.tangential_slip[component] =
                total_first * current_first[component] + total_second * current_second[component];
    }
    const NormalContactProperties& properties = _mechanical_properties[constraint.contact];
    result.pressure = constraint.friction_only ? equivalent_normal_pressure(constraint, state)
                                               : std::max(-properties.penalty * result.gap, 0.0);
    result.force = result.pressure * constraint.area;
    if (properties.friction_coefficient == 0.0 || !(result.pressure > 0.0)) return result;

    if (constraint.friction_only) {
        const std::array<double, 3> current_first = result.tangent_first,
                                    current_second = {normal[1] * current_first[2] - normal[2] * current_first[1],
                                        normal[2] * current_first[0] - normal[0] * current_first[2],
                                        normal[0] * current_first[1] - normal[1] * current_first[0]};
        std::array<double, 2> elastic_components{}, total_components{};
        for (std::size_t component = 0; component < 3; ++component) {
            elastic_components[0] += history.cartesian_elastic_tangential_slip[component] * current_first[component];
            elastic_components[1] += history.cartesian_elastic_tangential_slip[component] * current_second[component];
            total_components[0] += history.cartesian_total_tangential_slip[component] * current_first[component];
            total_components[1] += history.cartesian_total_tangential_slip[component] * current_second[component];
        }
        for (std::size_t node = 0; node < node_count; ++node)
            for (std::size_t component = 0; component < 3; ++component) {
                const double increment =
                    state[component * node_count + node] - committed_state[component * node_count + node];
                const std::size_t stored = stored_node(node);
                total_components[0] -= constraint.tangent_first_coefficients[stored][component] * increment;
                total_components[1] -= constraint.tangent_second_coefficients[stored][component] * increment;
                elastic_components[0] -= constraint.tangent_first_coefficients[stored][component] * increment;
                elastic_components[1] -= constraint.tangent_second_coefficients[stored][component] * increment;
            }
        for (std::size_t component = 0; component < 3; ++component) {
            result.tangential_slip[component] =
                total_components[0] * current_first[component] + total_components[1] * current_second[component];
            result.elastic_tangential_slip[component] =
                elastic_components[0] * current_first[component] + elastic_components[1] * current_second[component];
        }
        const std::array<double, 2> trial_elastic_components = elastic_components;
        const double sliding_limit = properties.friction_coefficient * result.pressure;
        result.stick_stiffness = properties.maximum_elastic_slip > 0.0 ? sliding_limit / properties.maximum_elastic_slip
                                                                       : properties.penalty;
        if (!std::isfinite(result.stick_stiffness) || !(result.stick_stiffness > 0.0))
            throw std::domain_error("Abaqus-style averaged contact has a nonpositive tangential stick stiffness");
        for (std::size_t tangent = 0; tangent < 2; ++tangent)
            result.trial_tangential_traction_components[tangent] = result.stick_stiffness * elastic_components[tangent];
        result.trial_tangential_magnitude =
            std::hypot(result.trial_tangential_traction_components[0], result.trial_tangential_traction_components[1]);
        if (result.trial_tangential_magnitude < sliding_limit ||
            (result.trial_tangential_magnitude == sliding_limit && !history.sliding)) {
            result.tangential_traction_components = result.trial_tangential_traction_components;
        } else {
            if (!(result.trial_tangential_magnitude > 0.0))
                throw std::domain_error("Abaqus-style averaged sliding contact has an undefined tangential direction");
            for (std::size_t tangent = 0; tangent < 2; ++tangent) {
                result.tangential_traction_components[tangent] = sliding_limit *
                                                                 result.trial_tangential_traction_components[tangent] /
                                                                 result.trial_tangential_magnitude;
                elastic_components[tangent] = result.tangential_traction_components[tangent] / result.stick_stiffness;
                result.friction_dissipation += constraint.area * result.tangential_traction_components[tangent] *
                                               (trial_elastic_components[tangent] - elastic_components[tangent]);
            }
            result.sliding = true;
            for (std::size_t component = 0; component < 3; ++component)
                result.elastic_tangential_slip[component] = elastic_components[0] * current_first[component] +
                                                            elastic_components[1] * current_second[component];
        }
        for (std::size_t component = 0; component < 3; ++component) {
            result.trial_tangential_traction[component] =
                result.trial_tangential_traction_components[0] * current_first[component] +
                result.trial_tangential_traction_components[1] * current_second[component];
            result.tangential_traction[component] =
                result.tangential_traction_components[0] * current_first[component] +
                result.tangential_traction_components[1] * current_second[component];
        }
        result.tangential_force = constraint.area * std::hypot(result.tangential_traction_components[0],
                                                        result.tangential_traction_components[1]);
        return result;
    }

    std::array<double, 3> relative_increment{};
    if (constraint.finite_sliding) {
        const std::array<double, 3> current_first = result.tangent_first,
                                    current_second = {normal[1] * current_first[2] - normal[2] * current_first[1],
                                        normal[2] * current_first[0] - normal[0] * current_first[2],
                                        normal[0] * current_first[1] - normal[1] * current_first[0]};
        std::array<double, 3> committed_normal = {constraint.reference_normal.x, constraint.reference_normal.y,
                                  constraint.reference_normal.z},
                              committed_first = {constraint.reference_tangent_first.x,
                                  constraint.reference_tangent_first.y, constraint.reference_tangent_first.z};
        if (history.cartesian_tangent_basis_initialized) {
            committed_normal = history.cartesian_contact_normal;
            committed_first = history.cartesian_contact_tangent_first;
        }
        const std::array<double, 3> committed_second = {
            committed_normal[1] * committed_first[2] - committed_normal[2] * committed_first[1],
            committed_normal[2] * committed_first[0] - committed_normal[0] * committed_first[2],
            committed_normal[0] * committed_first[1] - committed_normal[1] * committed_first[0]};
        std::array<double, 3> current_separation{}, committed_separation{}, transported_history{},
            transported_total_history{};
        for (std::size_t component = 0; component < 3; ++component) {
            for (std::size_t node = 0; node < node_count; ++node) {
                const std::size_t stored = stored_node(node);
                const double coordinate = component == 0   ? constraint.reference_coordinates[stored].x
                                          : component == 1 ? constraint.reference_coordinates[stored].y
                                                           : constraint.reference_coordinates[stored].z;
                current_separation[component] +=
                    constraint.gap_coefficients[stored] * (coordinate + state[component * node_count + node]);
                committed_separation[component] +=
                    constraint.gap_coefficients[stored] * (coordinate + committed_state[component * node_count + node]);
            }
        }
        double current_first_coordinate = 0.0, current_second_coordinate = 0.0, committed_first_coordinate = 0.0,
               committed_second_coordinate = 0.0, history_first = 0.0, history_second = 0.0;
        double total_history_first = 0.0, total_history_second = 0.0;
        for (std::size_t component = 0; component < 3; ++component) {
            current_first_coordinate -= current_separation[component] * current_first[component];
            current_second_coordinate -= current_separation[component] * current_second[component];
            committed_first_coordinate -= committed_separation[component] * committed_first[component];
            committed_second_coordinate -= committed_separation[component] * committed_second[component];
            history_first += history.cartesian_elastic_tangential_slip[component] * committed_first[component];
            history_second += history.cartesian_elastic_tangential_slip[component] * committed_second[component];
            total_history_first += history.cartesian_total_tangential_slip[component] * committed_first[component];
            total_history_second += history.cartesian_total_tangential_slip[component] * committed_second[component];
        }
        for (std::size_t component = 0; component < 3; ++component) {
            transported_history[component] =
                history_first * current_first[component] + history_second * current_second[component];
            transported_total_history[component] =
                total_history_first * current_first[component] + total_history_second * current_second[component];
            relative_increment[component] =
                (current_first_coordinate - committed_first_coordinate) * current_first[component] +
                (current_second_coordinate - committed_second_coordinate) * current_second[component];
            result.elastic_tangential_slip[component] = transported_history[component] + relative_increment[component];
            result.tangential_slip[component] = transported_total_history[component] + relative_increment[component];
        }
    } else {
        for (std::size_t component = 0; component < 3; ++component)
            relative_increment[component] = committed_relative[component] - relative[component];
    }
    double normal_increment = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        normal_increment += normal[component] * relative_increment[component];
    for (std::size_t component = 0; component < 3; ++component)
        if (!constraint.finite_sliding)
            result.elastic_tangential_slip[component] =
                history.cartesian_elastic_tangential_slip[component] + relative_increment[component];
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] -= normal_increment * normal[component];
    double history_normal = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        history_normal += result.elastic_tangential_slip[component] * normal[component];
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] -= history_normal * normal[component];
    const std::array<double, 3> trial_elastic_tangential_slip = result.elastic_tangential_slip;
    const double sliding_limit = properties.friction_coefficient * result.pressure;
    result.stick_stiffness =
        properties.maximum_elastic_slip > 0.0 ? sliding_limit / properties.maximum_elastic_slip : properties.penalty;
    if (!std::isfinite(result.stick_stiffness) || !(result.stick_stiffness > 0.0))
        throw std::domain_error("Abaqus-style averaged contact has a nonpositive tangential stick stiffness");
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
            throw std::domain_error("Abaqus-style averaged sliding contact has an undefined tangential direction");
        for (std::size_t component = 0; component < 3; ++component) {
            result.tangential_traction[component] =
                sliding_limit * result.trial_tangential_traction[component] / result.trial_tangential_magnitude;
            result.elastic_tangential_slip[component] = result.tangential_traction[component] / result.stick_stiffness;
        }
        result.sliding = true;
        for (std::size_t component = 0; component < 3; ++component)
            result.friction_dissipation +=
                constraint.area * result.tangential_traction[component] *
                (trial_elastic_tangential_slip[component] - result.elastic_tangential_slip[component]);
    }
    result.tangential_force = constraint.area * std::hypot(result.tangential_traction[0], result.tangential_traction[1],
                                                    result.tangential_traction[2]);
    return result;
}

SpatialAssembly::AbaqusAveragedConstraintValue SpatialAssembly::finite_region_normal_value(
    const AbaqusAveragedConstraint& constraint, const std::vector<std::size_t>& local_nodes,
    const std::vector<double>& state, std::vector<double>* residual, std::vector<double>* jacobian,
    std::vector<double>* pressure_derivative) const {
    if (!constraint.finite_region_normal)
        throw std::logic_error("HEX8 finite-region normal evaluation requires a finite-region constraint");
    const std::size_t node_count = local_nodes.size(), local_size = 3 * node_count;
    if (state.size() != local_size) throw std::invalid_argument("HEX8 finite-region normal state has the wrong size");
    const bool cached_state_matches =
        constraint.finite_region_cache_valid && constraint.finite_region_cached_state.size() == state.size() &&
        std::equal(state.begin(), state.end(), constraint.finite_region_cached_state.begin());
    if (residual == nullptr && jacobian == nullptr && cached_state_matches &&
        (pressure_derivative == nullptr || constraint.finite_region_cached_derivative_valid)) {
        if (pressure_derivative != nullptr) *pressure_derivative = constraint.finite_region_cached_pressure_derivative;
        AbaqusAveragedConstraintValue cached{};
        cached.gap = constraint.finite_region_cached_gap;
        cached.pressure = constraint.finite_region_cached_pressure;
        cached.force = constraint.finite_region_cached_force;
        return cached;
    }
    if (!constraint.projected) {
        if (residual != nullptr) residual->assign(local_size, 0.0);
        if (jacobian != nullptr) jacobian->assign(local_size * local_size, 0.0);
        if (pressure_derivative != nullptr) pressure_derivative->assign(local_size, 0.0);
        return {};
    }
    std::vector<double> unit_residual(local_size), area_derivative(local_size), gap_integral_derivative(local_size),
        unit_jacobian;
    if (jacobian != nullptr) unit_jacobian.resize(local_size * local_size);
    double area = 0.0, gap_integral = 0.0;
    for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
        const SecondaryContactFace& secondary =
            _secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
        for (const AbaqusAveragedConstraint::FiniteSlidingSample::NormalPoint& normal_point : sample.normal_points) {
            const PrimaryContactFace& primary =
                _primary_contact_faces.at(constraint.contact).at(normal_point.primary_face);
            const Quad4FaceQuadraturePoint point = make_quad4_face_quadrature_point(
                secondary.coordinates, normal_point.xi, normal_point.eta, normal_point.weight);
            const Quad4ToQuad4MechanicalGeometry geometry = {secondary.coordinates, primary.coordinates, point.shape,
                point.derivative_xi, point.derivative_eta, point.derivative_xi, point.derivative_eta,
                normal_point.weight, primary_material_orientation(primary.coordinates, primary.parent_centroid),
                sample.normal_orientation};
            Quad4SurfaceContactLocalValues local_state{};
            std::array<std::size_t, 8> local_node_indices{};
            for (std::size_t node = 0; node < 4; ++node) {
                const auto secondary_node = std::find(local_nodes.begin(), local_nodes.end(), secondary.nodes[node]);
                const auto primary_node = std::find(local_nodes.begin(), local_nodes.end(), primary.nodes[node]);
                if (secondary_node == local_nodes.end() || primary_node == local_nodes.end())
                    throw std::logic_error("HEX8 finite-region normal constraint lost a contact node");
                local_node_indices[node] = static_cast<std::size_t>(secondary_node - local_nodes.begin());
                local_node_indices[4 + node] = static_cast<std::size_t>(primary_node - local_nodes.begin());
            }
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t node = 0; node < 8; ++node)
                    local_state[8 * (component + 1) + node] = state[component * node_count + local_node_indices[node]];
            Quad4FiniteRegionNormalGeometryJacobian local_jacobian{};
            const Quad4FiniteRegionNormalGeometryValue value = compute_quad4_finite_region_normal_geometry(geometry,
                local_state, jacobian == nullptr && pressure_derivative == nullptr ? nullptr : &local_jacobian);
            if (!value.projected)
                throw std::domain_error("HEX8 finite-region normal integration point lost its active projection");
            area += value.area;
            gap_integral += value.gap_integral;
            for (std::size_t row_component = 0; row_component < 3; ++row_component)
                for (std::size_t row_node = 0; row_node < 8; ++row_node) {
                    const std::size_t local_row = 8 * (row_component + 1) + row_node,
                                      row = row_component * node_count + local_node_indices[row_node];
                    unit_residual[row] += value.unit_pressure_residual[local_row];
                    if (jacobian != nullptr)
                        for (std::size_t column_component = 0; column_component < 3; ++column_component)
                            for (std::size_t column_node = 0; column_node < 8; ++column_node) {
                                const std::size_t local_column = 8 * (column_component + 1) + column_node,
                                                  column =
                                                      column_component * node_count + local_node_indices[column_node];
                                unit_jacobian[row * local_size + column] +=
                                    local_jacobian[(2 + local_row) * quad4_surface_contact_local_dof_count +
                                                   local_column];
                            }
                }
            if (jacobian != nullptr || pressure_derivative != nullptr)
                for (std::size_t column_component = 0; column_component < 3; ++column_component)
                    for (std::size_t column_node = 0; column_node < 8; ++column_node) {
                        const std::size_t local_column = 8 * (column_component + 1) + column_node,
                                          column = column_component * node_count + local_node_indices[column_node];
                        area_derivative[column] += local_jacobian[local_column];
                        gap_integral_derivative[column] +=
                            local_jacobian[quad4_surface_contact_local_dof_count + local_column];
                    }
        }
    }
    if (!std::isfinite(area) || !(area > 0.0))
        throw std::domain_error("HEX8 finite-region normal constraint has a nonpositive current area");
    const double gap = gap_integral / area;
    const NormalContactProperties& properties = _mechanical_properties.at(constraint.contact);
    const double pressure = std::max(-properties.penalty * gap, 0.0);
    std::vector<double> local_pressure_derivative(local_size);
    if (pressure > 0.0)
        for (std::size_t column = 0; column < local_size; ++column)
            local_pressure_derivative[column] =
                -properties.penalty * (gap_integral_derivative[column] - gap * area_derivative[column]) / area;
    if (pressure_derivative != nullptr) *pressure_derivative = local_pressure_derivative;
    if (residual != nullptr) {
        residual->resize(local_size);
        for (std::size_t row = 0; row < local_size; ++row) (*residual)[row] = pressure * unit_residual[row];
    }
    if (jacobian != nullptr) {
        jacobian->resize(local_size * local_size);
        for (std::size_t row = 0; row < local_size; ++row)
            for (std::size_t column = 0; column < local_size; ++column)
                (*jacobian)[row * local_size + column] = pressure * unit_jacobian[row * local_size + column] +
                                                         unit_residual[row] * local_pressure_derivative[column];
    }
    constraint.finite_region_cached_state = state;
    constraint.finite_region_cached_pressure_derivative = local_pressure_derivative;
    constraint.finite_region_cached_gap = gap;
    constraint.finite_region_cached_pressure = pressure;
    constraint.finite_region_cached_force = pressure * area;
    constraint.finite_region_cache_valid = true;
    constraint.finite_region_cached_derivative_valid = jacobian != nullptr || pressure_derivative != nullptr;
    AbaqusAveragedConstraintValue result{};
    result.gap = gap;
    result.pressure = pressure;
    result.force = pressure * area;
    return result;
}

double SpatialAssembly::equivalent_normal_pressure(const AbaqusAveragedConstraint& constraint,
    const std::vector<double>& state, std::vector<double>* derivative,
    const std::vector<double>* friction_area_derivative) const {
    const auto finite_region = std::find_if(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
        [&constraint](const AbaqusAveragedConstraint& value) {
            return value.contact == constraint.contact && value.secondary == constraint.secondary &&
                   value.finite_region_normal;
        });
    if (finite_region != _abaqus_averaged_constraints.end()) {
        if (finite_region->nodes != constraint.nodes)
            throw std::logic_error("HEX8 averaged friction and finite-region normal supports differ");
        const std::vector<std::size_t>& local_nodes =
            constraint.active_nodes.empty() ? constraint.nodes : constraint.active_nodes;
        std::vector<double> pressure_derivative;
        const AbaqusAveragedConstraintValue value = finite_region_normal_value(*finite_region, local_nodes, state,
            nullptr, nullptr, derivative == nullptr ? nullptr : &pressure_derivative);
        if (derivative != nullptr) *derivative = std::move(pressure_derivative);
        return value.pressure;
    }
    double force = 0.0;
    std::vector<double> force_derivative(state.size());
    for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
        const auto point = std::find_if(
            _mechanical_points.begin(), _mechanical_points.end(), [&constraint, &sample](const MechanicalPoint& value) {
                return value.contact == constraint.contact && value.secondary_face == sample.secondary_face &&
                       value.secondary_local_point == sample.secondary_local_point;
            });
        if (point == _mechanical_points.end())
            throw std::logic_error("HEX8 averaged friction lost its matching normal constraint point");
        const std::size_t point_index = static_cast<std::size_t>(point - _mechanical_points.begin()),
                          primary = _mechanical_active_primary.at(point_index);
        if (primary == std::numeric_limits<std::size_t>::max())
            throw std::domain_error("HEX8 averaged friction lost its matching normal projection");
        const MechanicalCandidate candidate = mechanical_candidate(point_index, primary);
        if (!candidate.surface_to_surface)
            throw std::logic_error("HEX8 averaged friction requires surface-to-surface normal constraints");
        Quad4SurfaceContactLocalValues local_state{};
        for (std::size_t component = 0; component < 3; ++component)
            for (std::size_t node = 0; node < candidate.nodes.size(); ++node) {
                const std::vector<std::size_t>& local_nodes =
                    constraint.active_nodes.empty() ? constraint.nodes : constraint.active_nodes;
                const auto found = std::find(local_nodes.begin(), local_nodes.end(), candidate.nodes[node]);
                if (found == local_nodes.end())
                    throw std::logic_error("HEX8 averaged friction normal pressure lost a contact node");
                local_state[8 * (component + 1) + node] =
                    state[component * local_nodes.size() + static_cast<std::size_t>(found - local_nodes.begin())];
            }
        NormalContactProperties properties = _mechanical_properties.at(constraint.contact);
        properties.friction_coefficient = 0.0;
        properties.maximum_elastic_slip = 0.0;
        Quad4NormalForceAreaJacobian local_jacobian{};
        const Quad4NormalForceAreaValue value = compute_quad4_to_quad4_normal_force_area(
            properties, candidate.surface_geometry, local_state, derivative == nullptr ? nullptr : &local_jacobian);
        if (!value.projected) throw std::domain_error("HEX8 averaged friction lost its matching normal projection");
        force += value.force;
        if (derivative != nullptr)
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t node = 0; node < candidate.nodes.size(); ++node) {
                    const std::vector<std::size_t>& local_nodes =
                        constraint.active_nodes.empty() ? constraint.nodes : constraint.active_nodes;
                    const auto found = std::find(local_nodes.begin(), local_nodes.end(), candidate.nodes[node]);
                    const std::size_t local = 8 * (component + 1) + node,
                                      column = component * local_nodes.size() +
                                               static_cast<std::size_t>(found - local_nodes.begin());
                    force_derivative[column] += local_jacobian[local];
                }
    }
    if (!(constraint.area > 0.0)) {
        if (derivative != nullptr) derivative->assign(state.size(), 0.0);
        return 0.0;
    }
    const double pressure = force / constraint.area;
    if (derivative != nullptr) {
        if (friction_area_derivative == nullptr || friction_area_derivative->size() != state.size())
            throw std::logic_error("HEX8 averaged friction pressure requires its current-area derivative");
        derivative->resize(state.size());
        for (std::size_t column = 0; column < state.size(); ++column)
            (*derivative)[column] =
                (force_derivative[column] - pressure * friction_area_derivative->at(column)) / constraint.area;
    }
    return pressure;
}

void SpatialAssembly::compute_averaged_friction_geometry(const AbaqusAveragedConstraint& constraint,
    const std::vector<double>& state, const std::array<double, 2>& traction, std::vector<double>& residual,
    std::vector<double>* jacobian) const {
    const std::vector<std::size_t>& active_nodes =
        constraint.active_nodes.empty() ? constraint.nodes : constraint.active_nodes;
    const std::size_t node_count = active_nodes.size(), local_size = 3 * node_count;
    residual.assign(local_size, 0.0);
    if (jacobian != nullptr) jacobian->assign(local_size * local_size, 0.0);
    const Matrix4& averaging = abaqus_quad4_averaging();
    for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
        if (sample.primary_faces.empty()) continue;
        const SecondaryContactFace& secondary =
            _secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
        const std::array<double, 2>& location = abaqus_quad4_constraint_locations().at(sample.secondary_local_point);
        const Quad4FaceQuadraturePoint point = make_quad4_face_quadrature_point(
                                           secondary.coordinates, location[0], location[1], 1.0),
                                       normal_point = make_quad4_face_quadrature_point(secondary.coordinates,
                                           (4.0 / 3.0) * location[0], (4.0 / 3.0) * location[1], 1.0);
        std::array<double, 4> distribution{};
        for (std::size_t node = 0; node < 4; ++node)
            distribution[node] = averaging[sample.secondary_local_point * 4 + node];
        for (const std::size_t primary_index : sample.primary_faces) {
            const PrimaryContactFace& primary = _primary_contact_faces.at(constraint.contact).at(primary_index);
            const Quad4ToQuad4MechanicalGeometry geometry = {secondary.coordinates, primary.coordinates, point.shape,
                point.derivative_xi, point.derivative_eta, normal_point.derivative_xi, normal_point.derivative_eta,
                1.0 / static_cast<double>(sample.primary_faces.size()), 1.0, sample.normal_orientation};
            Quad4SurfaceContactLocalValues local_state{};
            std::array<std::size_t, 8> local_nodes{};
            for (std::size_t node = 0; node < 4; ++node) {
                const auto secondary_node = std::find(active_nodes.begin(), active_nodes.end(), secondary.nodes[node]);
                const auto primary_node = std::find(active_nodes.begin(), active_nodes.end(), primary.nodes[node]);
                if (secondary_node == active_nodes.end() || primary_node == active_nodes.end())
                    throw std::logic_error("HEX8 averaged friction geometry lost a contact node");
                local_nodes[node] = static_cast<std::size_t>(secondary_node - active_nodes.begin());
                local_nodes[4 + node] = static_cast<std::size_t>(primary_node - active_nodes.begin());
            }
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t node = 0; node < 8; ++node)
                    local_state[8 * (component + 1) + node] = state[component * node_count + local_nodes[node]];
            Quad4SurfaceContactLocalJacobian local_jacobian{};
            const Quad4SurfaceContactLocalResidual local_residual =
                compute_quad4_to_quad4_tangential_force_geometry(geometry, distribution, sample.tangent_orientation,
                    traction[0], traction[1], local_state, jacobian == nullptr ? nullptr : &local_jacobian);
            for (std::size_t row_component = 0; row_component < 3; ++row_component)
                for (std::size_t row_node = 0; row_node < 8; ++row_node) {
                    const std::size_t local_row = 8 * (row_component + 1) + row_node,
                                      row = row_component * node_count + local_nodes[row_node];
                    residual[row] += local_residual[local_row];
                    if (jacobian != nullptr)
                        for (std::size_t column_component = 0; column_component < 3; ++column_component)
                            for (std::size_t column_node = 0; column_node < 8; ++column_node) {
                                const std::size_t local_column = 8 * (column_component + 1) + column_node,
                                                  column = column_component * node_count + local_nodes[column_node];
                                (*jacobian)[row * local_size + column] +=
                                    local_jacobian[local_row * quad4_surface_contact_local_dof_count + local_column];
                            }
                }
        }
    }
}

void SpatialAssembly::compute_averaged_friction_traction_derivatives(const AbaqusAveragedConstraint& constraint,
    const std::vector<double>& state, const std::vector<double>& committed_state, const ContactPointHistory& history,
    const AbaqusAveragedConstraintValue& value, std::array<std::vector<double>, 2>& derivatives) const {
    const std::vector<std::size_t>& active_nodes =
        constraint.active_nodes.empty() ? constraint.nodes : constraint.active_nodes;
    const std::size_t node_count = active_nodes.size(), local_size = 3 * node_count;
    double area = 0.0;
    std::array<double, 3> normal_raw{}, tangent_raw{};
    std::array<double, 2> increment{};
    std::vector<double> area_derivative(local_size);
    std::array<std::vector<double>, 3> normal_raw_derivative, tangent_raw_derivative;
    std::array<std::vector<double>, 2> increment_derivative;
    for (std::size_t component = 0; component < 3; ++component) {
        normal_raw_derivative[component].resize(local_size);
        tangent_raw_derivative[component].resize(local_size);
    }
    for (std::size_t tangent = 0; tangent < 2; ++tangent) increment_derivative[tangent].resize(local_size);
    const Matrix4& averaging = abaqus_quad4_averaging();
    for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
        if (sample.primary_faces.empty()) continue;
        const SecondaryContactFace& secondary =
            _secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
        const std::array<double, 2>& location = abaqus_quad4_constraint_locations().at(sample.secondary_local_point);
        const Quad4FaceQuadraturePoint point = make_quad4_face_quadrature_point(
                                           secondary.coordinates, location[0], location[1], 1.0),
                                       normal_point = make_quad4_face_quadrature_point(secondary.coordinates,
                                           (4.0 / 3.0) * location[0], (4.0 / 3.0) * location[1], 1.0);
        std::array<double, 4> distribution{};
        for (std::size_t node = 0; node < 4; ++node)
            distribution[node] = averaging[sample.secondary_local_point * 4 + node];
        for (const std::size_t primary_index : sample.primary_faces) {
            const PrimaryContactFace& primary = _primary_contact_faces.at(constraint.contact).at(primary_index);
            const Quad4ToQuad4MechanicalGeometry geometry = {secondary.coordinates, primary.coordinates, point.shape,
                point.derivative_xi, point.derivative_eta, normal_point.derivative_xi, normal_point.derivative_eta,
                1.0 / static_cast<double>(sample.primary_faces.size()), 1.0, sample.normal_orientation};
            Quad4SurfaceContactLocalValues local_state{}, local_committed{};
            std::array<std::size_t, 8> local_nodes{};
            for (std::size_t node = 0; node < 4; ++node) {
                const auto secondary_node = std::find(active_nodes.begin(), active_nodes.end(), secondary.nodes[node]);
                const auto primary_node = std::find(active_nodes.begin(), active_nodes.end(), primary.nodes[node]);
                if (secondary_node == active_nodes.end() || primary_node == active_nodes.end())
                    throw std::logic_error("HEX8 averaged friction derivative lost a contact node");
                local_nodes[node] = static_cast<std::size_t>(secondary_node - active_nodes.begin());
                local_nodes[4 + node] = static_cast<std::size_t>(primary_node - active_nodes.begin());
            }
            for (std::size_t component = 0; component < 3; ++component)
                for (std::size_t node = 0; node < 8; ++node) {
                    const std::size_t source = component * node_count + local_nodes[node],
                                      target = 8 * (component + 1) + node;
                    local_state[target] = state[source];
                    local_committed[target] = committed_state[source];
                }
            Quad4AveragedFrictionGeometryJacobian local_jacobian{};
            const Quad4AveragedFrictionGeometryValue geometry_value = compute_quad4_averaged_friction_geometry_value(
                geometry, distribution, sample.tangent_orientation, local_state, local_committed, &local_jacobian);
            if (!geometry_value.projected)
                throw std::domain_error("HEX8 averaged friction derivative lost its selected primary face");
            area += geometry_value.area;
            for (std::size_t component = 0; component < 3; ++component) {
                normal_raw[component] += geometry_value.area_normal[component];
                tangent_raw[component] += geometry_value.area_tangent_first[component];
            }
            for (std::size_t tangent = 0; tangent < 2; ++tangent)
                increment[tangent] += geometry_value.tangential_increment[tangent];
            for (std::size_t column_component = 0; column_component < 3; ++column_component)
                for (std::size_t column_node = 0; column_node < 8; ++column_node) {
                    const std::size_t local_column = 8 * (column_component + 1) + column_node,
                                      column = column_component * node_count + local_nodes[column_node];
                    area_derivative[column] += local_jacobian[local_column];
                    for (std::size_t component = 0; component < 3; ++component) {
                        normal_raw_derivative[component][column] +=
                            local_jacobian[(1 + component) * quad4_surface_contact_local_dof_count + local_column];
                        tangent_raw_derivative[component][column] +=
                            local_jacobian[(4 + component) * quad4_surface_contact_local_dof_count + local_column];
                    }
                    for (std::size_t tangent = 0; tangent < 2; ++tangent)
                        increment_derivative[tangent][column] +=
                            local_jacobian[(10 + tangent) * quad4_surface_contact_local_dof_count + local_column];
                }
        }
    }
    if (!(area > 0.0)) throw std::domain_error("HEX8 averaged friction derivative has a nonpositive area");
    const double normal_measure =
        std::sqrt(normal_raw[0] * normal_raw[0] + normal_raw[1] * normal_raw[1] + normal_raw[2] * normal_raw[2]);
    std::array<double, 3> normal{};
    std::array<std::vector<double>, 3> normal_derivative;
    for (std::size_t component = 0; component < 3; ++component) {
        normal[component] = normal_raw[component] / normal_measure;
        normal_derivative[component].resize(local_size);
    }
    for (std::size_t column = 0; column < local_size; ++column) {
        double measure_derivative = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            measure_derivative += normal[component] * normal_raw_derivative[component][column];
        for (std::size_t component = 0; component < 3; ++component)
            normal_derivative[component][column] =
                (normal_raw_derivative[component][column] - normal[component] * measure_derivative) / normal_measure;
    }
    double tangent_normal = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        tangent_normal += tangent_raw[component] * normal[component];
    std::array<double, 3> tangent_projected{};
    std::array<std::vector<double>, 3> tangent_projected_derivative;
    for (std::size_t component = 0; component < 3; ++component) {
        tangent_projected[component] = tangent_raw[component] - tangent_normal * normal[component];
        tangent_projected_derivative[component].resize(local_size);
    }
    for (std::size_t column = 0; column < local_size; ++column) {
        double tangent_normal_derivative = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            tangent_normal_derivative += tangent_raw_derivative[component][column] * normal[component] +
                                         tangent_raw[component] * normal_derivative[component][column];
        for (std::size_t component = 0; component < 3; ++component)
            tangent_projected_derivative[component][column] = tangent_raw_derivative[component][column] -
                                                              tangent_normal_derivative * normal[component] -
                                                              tangent_normal * normal_derivative[component][column];
    }
    const double tangent_measure =
        std::sqrt(tangent_projected[0] * tangent_projected[0] + tangent_projected[1] * tangent_projected[1] +
                  tangent_projected[2] * tangent_projected[2]);
    std::array<double, 3> first{}, second{};
    std::array<std::vector<double>, 3> first_derivative, second_derivative;
    for (std::size_t component = 0; component < 3; ++component) {
        first[component] = tangent_projected[component] / tangent_measure;
        first_derivative[component].resize(local_size);
        second_derivative[component].resize(local_size);
    }
    second = {normal[1] * first[2] - normal[2] * first[1], normal[2] * first[0] - normal[0] * first[2],
        normal[0] * first[1] - normal[1] * first[0]};
    for (std::size_t column = 0; column < local_size; ++column) {
        double tangent_measure_derivative = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            tangent_measure_derivative += first[component] * tangent_projected_derivative[component][column];
        for (std::size_t component = 0; component < 3; ++component)
            first_derivative[component][column] =
                (tangent_projected_derivative[component][column] - first[component] * tangent_measure_derivative) /
                tangent_measure;
        second_derivative[0][column] =
            normal_derivative[1][column] * first[2] + normal[1] * first_derivative[2][column] -
            normal_derivative[2][column] * first[1] - normal[2] * first_derivative[1][column];
        second_derivative[1][column] =
            normal_derivative[2][column] * first[0] + normal[2] * first_derivative[0][column] -
            normal_derivative[0][column] * first[2] - normal[0] * first_derivative[2][column];
        second_derivative[2][column] =
            normal_derivative[0][column] * first[1] + normal[0] * first_derivative[1][column] -
            normal_derivative[1][column] * first[0] - normal[1] * first_derivative[0][column];
    }
    std::array<double, 2> elastic{};
    std::array<std::vector<double>, 2> elastic_derivative;
    const std::array<double, 3>& stored = history.cartesian_elastic_tangential_slip;
    for (std::size_t tangent = 0; tangent < 2; ++tangent) {
        elastic_derivative[tangent].resize(local_size);
        const std::array<double, 3>& basis = tangent == 0 ? first : second;
        const std::array<std::vector<double>, 3>& basis_derivative =
            tangent == 0 ? first_derivative : second_derivative;
        for (std::size_t component = 0; component < 3; ++component)
            elastic[tangent] += stored[component] * basis[component];
        elastic[tangent] -= increment[tangent] / area;
        for (std::size_t column = 0; column < local_size; ++column) {
            for (std::size_t component = 0; component < 3; ++component)
                elastic_derivative[tangent][column] += stored[component] * basis_derivative[component][column];
            elastic_derivative[tangent][column] -=
                (increment_derivative[tangent][column] - increment[tangent] / area * area_derivative[column]) / area;
        }
    }
    const NormalContactProperties& properties = _mechanical_properties[constraint.contact];
    std::vector<double> pressure_derivative;
    const double pressure = equivalent_normal_pressure(constraint, state, &pressure_derivative, &area_derivative),
                 sliding_limit = properties.friction_coefficient * pressure,
                 stick_stiffness = properties.maximum_elastic_slip > 0.0
                                       ? sliding_limit / properties.maximum_elastic_slip
                                       : properties.penalty;
    std::vector<double> stick_derivative(local_size);
    for (std::size_t column = 0; column < local_size; ++column)
        if (properties.maximum_elastic_slip > 0.0)
            stick_derivative[column] =
                properties.friction_coefficient * pressure_derivative[column] / properties.maximum_elastic_slip;
    std::array<double, 2> trial{};
    std::array<std::vector<double>, 2> trial_derivative;
    for (std::size_t tangent = 0; tangent < 2; ++tangent) {
        trial[tangent] = stick_stiffness * elastic[tangent];
        trial_derivative[tangent].resize(local_size);
        for (std::size_t column = 0; column < local_size; ++column)
            trial_derivative[tangent][column] =
                stick_stiffness * elastic_derivative[tangent][column] + stick_derivative[column] * elastic[tangent];
    }
    const double trial_magnitude = std::hypot(trial[0], trial[1]);
    derivatives[0].assign(local_size, 0.0);
    derivatives[1].assign(local_size, 0.0);
    if (!value.sliding) {
        derivatives = trial_derivative;
        return;
    }
    for (std::size_t column = 0; column < local_size; ++column) {
        const double magnitude_derivative =
                         (trial[0] * trial_derivative[0][column] + trial[1] * trial_derivative[1][column]) /
                         trial_magnitude,
                     limit_derivative = properties.friction_coefficient * pressure_derivative[column];
        for (std::size_t tangent = 0; tangent < 2; ++tangent)
            derivatives[tangent][column] =
                limit_derivative * trial[tangent] / trial_magnitude +
                sliding_limit * (trial_derivative[tangent][column] / trial_magnitude -
                                    trial[tangent] * magnitude_derivative / (trial_magnitude * trial_magnitude));
    }
}

void SpatialAssembly::compute_averaged_constraint(const AbaqusAveragedConstraint& constraint,
    const std::vector<double>& state, std::vector<double>& residual, std::vector<double>* jacobian) const {
    if (constraint.finite_region_normal) {
        finite_region_normal_value(constraint, constraint.active_nodes, state, &residual, jacobian);
        return;
    }
    const std::size_t node_count =
                          constraint.active_nodes.empty() ? constraint.nodes.size() : constraint.active_nodes.size(),
                      local_size = 3 * node_count;
    const auto stored_node = [&constraint](std::size_t local) {
        return constraint.active_node_indices.empty() ? local : constraint.active_node_indices.at(local);
    };
    std::vector<std::size_t> dofs;
    averaged_constraint_dofs(constraint, dofs);
    std::vector<double> committed_state(local_size);
    for (std::size_t local = 0; local < local_size; ++local)
        committed_state[local] = _committed_contact_solution.at(dofs[local]);
    const ContactPointHistory& history = _contact_histories[constraint.contact][constraint.history];
    const AbaqusAveragedConstraintValue value = averaged_constraint_value(constraint, state, committed_state, history);
    residual.assign(local_size, 0.0);
    if (jacobian != nullptr) jacobian->assign(local_size * local_size, 0.0);
    if (!(value.pressure > 0.0)) return;
    const std::array<double, 3> normal = {constraint.normal.x, constraint.normal.y, constraint.normal.z};
    const bool local_normal_gradients = !constraint.normal_gap_coefficients.empty();
    const auto gap_gradient = [&constraint, &normal, &stored_node, local_normal_gradients](
                                  std::size_t node, std::size_t component) {
        const std::size_t stored = stored_node(node);
        return local_normal_gradients ? constraint.normal_gap_coefficients[stored][component]
                                      : constraint.gap_coefficients[stored] * normal[component];
    };
    if (constraint.friction_only) {
        for (std::size_t row_component = 0; row_component < 3; ++row_component)
            for (std::size_t row_node = 0; row_node < node_count; ++row_node) {
                const std::size_t row = row_component * node_count + row_node, stored = stored_node(row_node);
                residual[row] = -constraint.area * (constraint.traction_first_coefficients[stored][row_component] *
                                                           value.tangential_traction_components[0] +
                                                       constraint.traction_second_coefficients[stored][row_component] *
                                                           value.tangential_traction_components[1]);
            }
        if (jacobian == nullptr) return;
        std::vector<double> geometric_residual;
        compute_averaged_friction_geometry(
            constraint, state, value.tangential_traction_components, geometric_residual, jacobian);
        std::array<std::vector<double>, 2> traction_derivatives;
        compute_averaged_friction_traction_derivatives(
            constraint, state, committed_state, history, value, traction_derivatives);
        for (std::size_t row_component = 0; row_component < 3; ++row_component)
            for (std::size_t row_node = 0; row_node < node_count; ++row_node) {
                const std::size_t row = row_component * node_count + row_node;
                const std::size_t stored = stored_node(row_node);
                const std::array<double, 2> row_coefficient = {
                    constraint.traction_first_coefficients[stored][row_component],
                    constraint.traction_second_coefficients[stored][row_component]};
                for (std::size_t column = 0; column < local_size; ++column)
                    for (std::size_t tangent = 0; tangent < 2; ++tangent)
                        (*jacobian)[row * local_size + column] -=
                            constraint.area * row_coefficient[tangent] * traction_derivatives[tangent][column];
            }
        return;
    }
    for (std::size_t component = 0; component < 3; ++component) {
        for (std::size_t node = 0; node < node_count; ++node) {
            const std::size_t local = component * node_count + node;
            residual[local] = -constraint.area * (value.pressure * gap_gradient(node, component) +
                                                     constraint.gap_coefficients[stored_node(node)] *
                                                         value.tangential_traction[component]);
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
            const double row_coefficient = constraint.gap_coefficients[stored_node(row_node)];
            for (std::size_t column_component = 0; column_component < 3; ++column_component)
                for (std::size_t column_node = 0; column_node < node_count; ++column_node) {
                    const std::size_t column = column_component * node_count + column_node;
                    const double column_coefficient = constraint.gap_coefficients[stored_node(column_node)],
                                 tangential_column_coefficient = -column_coefficient,
                                 pressure_derivative =
                                     -properties.penalty * gap_gradient(column_node, column_component);
                    double traction_derivative = local_normal_gradients || constraint.friction_only
                                                     ? 0.0
                                                     : normal[row_component] * pressure_derivative;
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
                    (*jacobian)[row * local_size + column] =
                        -constraint.area *
                        ((local_normal_gradients ? gap_gradient(row_node, row_component) * pressure_derivative : 0.0) +
                            row_coefficient * traction_derivative);
                }
        }
}

void SpatialAssembly::refresh_hex20_finite_averaged_constraints(const std::vector<double>& state) const {
    const auto current_face = [this, &state](
                                  const std::array<std::size_t, 8>& nodes, const Quad8FaceCoordinates& reference) {
        Quad8FaceCoordinates result = reference;
        for (std::size_t node = 0; node < 8; ++node) {
            const double displacement_x = state.at(dof(Field::displacement_x, nodes[node]));
            const double displacement_y = state.at(dof(Field::displacement_y, nodes[node]));
            const double displacement_z = state.at(dof(Field::displacement_z, nodes[node]));
            if (!std::isfinite(displacement_x) || !std::isfinite(displacement_y) || !std::isfinite(displacement_z))
                throw std::logic_error("HEX20 finite-sliding averaged constraint is missing the required shadow "
                                       "displacement of global node " +
                                       std::to_string(nodes[node]));
            result[node].x += displacement_x;
            result[node].y += displacement_y;
            result[node].z += displacement_z;
        }
        return result;
    };
    const Matrix8& averaging = abaqus_quad8_averaging();
    for (AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        if (!constraint.finite_sliding) continue;
        if (_touched_mechanical_nodes[mechanical_node_index(constraint.contact, constraint.history)] == 0U) continue;
        std::vector<std::size_t> required_nodes;
        std::fill(constraint.gap_coefficients.begin(), constraint.gap_coefficients.end(), 0.0);
        std::fill(constraint.secondary_coefficients.begin(), constraint.secondary_coefficients.end(), 0.0);
        std::fill(constraint.normal_gap_coefficients.begin(), constraint.normal_gap_coefficients.end(),
            std::array<double, 3>{});
        std::fill(constraint.secondary_normal_coefficients.begin(), constraint.secondary_normal_coefficients.end(),
            std::array<double, 3>{});
        constraint.normal = {};
        constraint.tangent_first = {};
        constraint.area = 0.0;
        constraint.projected = true;
        std::vector<double> primary_face_weight(_hex20_primary_contact_faces.at(constraint.contact).size(), 0.0);
        for (AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
            sample.primary_faces.clear();
            const Hex20SecondaryContactFace& secondary =
                _hex20_secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
            required_nodes.insert(
                required_nodes.end(), secondary.displacement_nodes.begin(), secondary.displacement_nodes.end());
            const Quad8FaceCoordinates secondary_current =
                current_face(secondary.displacement_nodes, secondary.coordinates);
            const Quad8FaceGeometry current_geometry = make_quad8_face_geometry(secondary_current);
            double face_area = 0.0;
            for (const Quad8FaceMechanicalQuadraturePoint& point : current_geometry.mechanical_points)
                face_area += point.quadrature_weight * reference_measure(point);
            const double fraction = sample.secondary_local_point < 4 ? 1.0 / 24.0 : 5.0 / 24.0,
                         local_area = face_area * fraction;
            if (!std::isfinite(local_area) || !(local_area > 0.0))
                throw std::domain_error("HEX20 finite-sliding averaged constraint has a nonpositive current area");
            constraint.area += local_area;
            const std::array<double, 2>& location =
                abaqus_quad8_constraint_locations().at(sample.secondary_local_point);
            const Quad8FaceMechanicalQuadraturePoint constraint_point =
                make_quad8_face_mechanical_point(secondary_current, location[0], location[1], 1.0);
            const CartesianPoint3 area_vector = cross(constraint_point.tangent_xi, constraint_point.tangent_eta);
            const double normal_measure = std::sqrt(dot(area_vector, area_vector));
            if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
                throw std::domain_error("HEX20 finite-sliding averaged constraint has an undefined current normal");
            const CartesianPoint3 local_normal = {sample.normal_orientation * area_vector.x / normal_measure,
                sample.normal_orientation * area_vector.y / normal_measure,
                sample.normal_orientation * area_vector.z / normal_measure};
            for (std::size_t local_node = 0; local_node < 8; ++local_node) {
                const double coefficient = local_area * averaging[sample.secondary_local_point * 8 + local_node];
                const auto node = std::find(
                    constraint.nodes.begin(), constraint.nodes.end(), secondary.displacement_nodes[local_node]);
                if (node == constraint.nodes.end())
                    throw std::logic_error("HEX20 finite-sliding averaged constraint lost a secondary node");
                const std::size_t stored_node = static_cast<std::size_t>(node - constraint.nodes.begin());
                constraint.gap_coefficients[stored_node] -= coefficient;
                constraint.normal_gap_coefficients[stored_node][0] -= coefficient * local_normal.x;
                constraint.normal_gap_coefficients[stored_node][1] -= coefficient * local_normal.y;
                constraint.normal_gap_coefficients[stored_node][2] -= coefficient * local_normal.z;
                const auto output =
                    std::find(_hex20_secondary_boundaries.at(constraint.contact).boundary.displacement_nodes.begin(),
                        _hex20_secondary_boundaries.at(constraint.contact).boundary.displacement_nodes.end(),
                        _hex20_secondary_boundaries.at(constraint.contact)
                            .boundary.faces.at(sample.secondary_face)
                            .nodes[local_node]);
                if (output == _hex20_secondary_boundaries.at(constraint.contact).boundary.displacement_nodes.end())
                    throw std::logic_error("HEX20 finite-sliding averaged constraint lost an output node");
                const std::size_t output_index = static_cast<std::size_t>(
                    output - _hex20_secondary_boundaries.at(constraint.contact).boundary.displacement_nodes.begin());
                const auto stored = std::find(
                    constraint.secondary_output_nodes.begin(), constraint.secondary_output_nodes.end(), output_index);
                if (stored == constraint.secondary_output_nodes.end())
                    throw std::logic_error("HEX20 finite-sliding averaged constraint output support changed");
                const std::size_t stored_output =
                    static_cast<std::size_t>(stored - constraint.secondary_output_nodes.begin());
                constraint.secondary_coefficients[stored_output] += coefficient;
                constraint.secondary_normal_coefficients[stored_output][0] += coefficient * local_normal.x;
                constraint.secondary_normal_coefficients[stored_output][1] += coefficient * local_normal.y;
                constraint.secondary_normal_coefficients[stored_output][2] += coefficient * local_normal.z;
            }
            constraint.normal.x += local_area * local_normal.x;
            constraint.normal.y += local_area * local_normal.y;
            constraint.normal.z += local_area * local_normal.z;
            CartesianPoint3 tangent = constraint_point.tangent_xi;
            const double tangent_measure = std::sqrt(dot(tangent, tangent));
            if (!std::isfinite(tangent_measure) || !(tangent_measure > 0.0))
                throw std::domain_error("HEX20 finite-sliding averaged constraint has an undefined current tangent");
            constraint.tangent_first.x += local_area * sample.tangent_orientation * tangent.x / tangent_measure;
            constraint.tangent_first.y += local_area * sample.tangent_orientation * tangent.y / tangent_measure;
            constraint.tangent_first.z += local_area * sample.tangent_orientation * tangent.z / tangent_measure;

            const auto& primary_faces = _hex20_primary_contact_faces.at(constraint.contact);
            for (const AbaqusQuad8TransferSample& transfer :
                abaqus_quad8_primary_transfer_rule(sample.secondary_local_point)) {
                const double xi = 2.0 * transfer.first - 1.0, eta = 2.0 * transfer.second - 1.0;
                const Quad8FaceMechanicalQuadraturePoint secondary_point =
                    make_quad8_face_mechanical_point(secondary_current, xi, eta, 1.0);
                double minimum_distance = std::numeric_limits<double>::infinity();
                std::size_t selected = std::numeric_limits<std::size_t>::max();
                Quad8ReferenceProjectionValue selected_projection{};
                const auto consider = [&](std::size_t primary_index) {
                    const Hex20PrimaryContactFace& primary = primary_faces[primary_index];
                    const Quad8FaceCoordinates primary_current =
                        current_face(primary.displacement_nodes, primary.coordinates);
                    const Quad8ReferenceProjectionValue projection = compute_quad8_reference_projection(
                        secondary_current, primary_current, secondary_point.displacement_shape,
                        normal_orientation(primary.coordinates, secondary.parent_centroid, primary.parent_centroid));
                    if (!projection.projected) return;
                    const double distance = std::abs(projection.gap);
                    if (distance < minimum_distance || (distance == minimum_distance && primary_index < selected)) {
                        minimum_distance = distance;
                        selected = primary_index;
                        selected_projection = projection;
                    }
                };
                if (primary_faces.size() > spatial_detail::contact_search_tree_minimum_items) {
                    std::array<double, 3> search_point{};
                    for (std::size_t node = 0; node < 8; ++node) {
                        search_point[0] += secondary_point.displacement_shape[node] * secondary_current[node].x;
                        search_point[1] += secondary_point.displacement_shape[node] * secondary_current[node].y;
                        search_point[2] += secondary_point.displacement_shape[node] * secondary_current[node].z;
                    }
                    _contact_search_trees[constraint.contact].begin_query(search_point, _contact_search_query);
                    std::size_t primary_index = 0;
                    while (_contact_search_trees[constraint.contact].next_candidate(
                        _contact_search_query, minimum_distance, primary_index))
                        consider(primary_index);
                } else {
                    for (std::size_t primary_index = 0; primary_index < primary_faces.size(); ++primary_index)
                        consider(primary_index);
                }
                if (selected == std::numeric_limits<std::size_t>::max()) {
                    constraint.projected = false;
                    continue;
                }
                if (std::find(sample.primary_faces.begin(), sample.primary_faces.end(), selected) ==
                    sample.primary_faces.end())
                    sample.primary_faces.push_back(selected);
                primary_face_weight[selected] += local_area * transfer.weight;
                const Hex20PrimaryContactFace& primary = primary_faces[selected];
                required_nodes.insert(
                    required_nodes.end(), primary.displacement_nodes.begin(), primary.displacement_nodes.end());
                for (std::size_t local_node = 0; local_node < 8; ++local_node) {
                    const auto node = std::find(
                        constraint.nodes.begin(), constraint.nodes.end(), primary.displacement_nodes[local_node]);
                    if (node == constraint.nodes.end())
                        throw std::logic_error("HEX20 finite-sliding averaged constraint lost a primary node");
                    const std::size_t stored_node = static_cast<std::size_t>(node - constraint.nodes.begin());
                    const double coefficient =
                        local_area * transfer.weight * selected_projection.primary_shape[local_node];
                    constraint.gap_coefficients[stored_node] += coefficient;
                    constraint.normal_gap_coefficients[stored_node][0] += coefficient * local_normal.x;
                    constraint.normal_gap_coefficients[stored_node][1] += coefficient * local_normal.y;
                    constraint.normal_gap_coefficients[stored_node][2] += coefficient * local_normal.z;
                }
            }
        }
        if (!std::isfinite(constraint.area) || !(constraint.area > 0.0))
            throw std::domain_error("HEX20 finite-sliding averaged constraint has a nonpositive current area");
        const double normal_measure = std::sqrt(dot(constraint.normal, constraint.normal));
        if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
            throw std::domain_error(
                "HEX20 finite-sliding averaged constraint has an undefined averaged current normal");
        constraint.normal.x /= normal_measure;
        constraint.normal.y /= normal_measure;
        constraint.normal.z /= normal_measure;
        const double tangent_normal_component = dot(constraint.tangent_first, constraint.normal);
        constraint.tangent_first.x -= tangent_normal_component * constraint.normal.x;
        constraint.tangent_first.y -= tangent_normal_component * constraint.normal.y;
        constraint.tangent_first.z -= tangent_normal_component * constraint.normal.z;
        const double tangent_measure = std::sqrt(dot(constraint.tangent_first, constraint.tangent_first));
        if (!std::isfinite(tangent_measure) || !(tangent_measure > 0.0))
            throw std::domain_error(
                "HEX20 finite-sliding averaged constraint has an undefined averaged current tangent");
        constraint.tangent_first.x /= tangent_measure;
        constraint.tangent_first.y /= tangent_measure;
        constraint.tangent_first.z /= tangent_measure;
        for (double& coefficient : constraint.gap_coefficients) coefficient /= constraint.area;
        for (double& coefficient : constraint.secondary_coefficients) coefficient /= constraint.area;
        for (std::array<double, 3>& coefficient : constraint.normal_gap_coefficients)
            for (double& component : coefficient) component /= constraint.area;
        for (std::array<double, 3>& coefficient : constraint.secondary_normal_coefficients)
            for (double& component : coefficient) component /= constraint.area;
        if (!primary_face_weight.empty())
            constraint.primary_face = static_cast<std::size_t>(
                std::max_element(primary_face_weight.begin(), primary_face_weight.end()) - primary_face_weight.begin());
        constraint.reference_gap = 0.0;
        for (std::size_t node = 0; node < constraint.nodes.size(); ++node) {
            constraint.reference_gap +=
                constraint.normal_gap_coefficients[node][0] * constraint.reference_coordinates[node].x +
                constraint.normal_gap_coefficients[node][1] * constraint.reference_coordinates[node].y +
                constraint.normal_gap_coefficients[node][2] * constraint.reference_coordinates[node].z;
        }
        std::sort(required_nodes.begin(), required_nodes.end());
        required_nodes.erase(std::unique(required_nodes.begin(), required_nodes.end()), required_nodes.end());
        constraint.active_nodes.clear();
        constraint.active_node_indices.clear();
        for (std::size_t stored = 0; stored < constraint.nodes.size(); ++stored)
            if (std::binary_search(required_nodes.begin(), required_nodes.end(), constraint.nodes[stored])) {
                constraint.active_nodes.push_back(constraint.nodes[stored]);
                constraint.active_node_indices.push_back(stored);
            }
        if (constraint.active_nodes.empty())
            throw std::logic_error("HEX20 finite-sliding averaged constraint has no active nodes");
    }
}

void SpatialAssembly::refresh_finite_averaged_constraints(const std::vector<double>& state) const {
    if (_uses_hex20) {
        refresh_hex20_finite_averaged_constraints(state);
        return;
    }
    const auto current_face = [this, &state](
                                  const std::array<std::size_t, 4>& nodes, const Quad4FaceCoordinates& reference) {
        Quad4FaceCoordinates result = reference;
        for (std::size_t node = 0; node < 4; ++node) {
            result[node].x += state.at(dof(Field::displacement_x, nodes[node]));
            result[node].y += state.at(dof(Field::displacement_y, nodes[node]));
            result[node].z += state.at(dof(Field::displacement_z, nodes[node]));
        }
        return result;
    };
    std::vector<std::vector<Quad4FaceCoordinates>> primary_current_faces_by_contact(_definition.contacts.size());
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        const bool averaged = std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
            [contact](const AbaqusAveragedConstraint& constraint) {
                return constraint.contact == contact && constraint.finite_sliding;
            });
        if (!averaged) continue;
        const auto& primary_faces = _primary_contact_faces.at(contact);
        auto& primary_current_faces = primary_current_faces_by_contact[contact];
        primary_current_faces.reserve(primary_faces.size());
        for (const PrimaryContactFace& face : primary_faces)
            primary_current_faces.push_back(current_face(face.nodes, face.coordinates));
    }
    const Matrix4& averaging = abaqus_quad4_averaging();
    constexpr double boundary_tolerance = 1.0e-9;
    for (AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        constraint.finite_region_cache_valid = false;
        constraint.finite_region_cached_derivative_valid = false;
        if (!constraint.finite_sliding) continue;
        if (constraint.friction_only) {
            const auto finite_region = std::find_if(_abaqus_averaged_constraints.begin(),
                _abaqus_averaged_constraints.end(), [&constraint](const AbaqusAveragedConstraint& value) {
                    return value.contact == constraint.contact && value.secondary == constraint.secondary &&
                           value.finite_region_normal;
                });
            if (finite_region != _abaqus_averaged_constraints.end()) {
                constraint.finite_sliding_samples = finite_region->finite_sliding_samples;
                constraint.normal = finite_region->normal;
                constraint.tangent_first = finite_region->tangent_first;
                constraint.projected = finite_region->projected;
                constraint.gap_coefficients = finite_region->gap_coefficients;
                constraint.secondary_coefficients = finite_region->secondary_coefficients;
                constraint.tangent_first_coefficients = finite_region->tangent_first_coefficients;
                constraint.tangent_second_coefficients = finite_region->tangent_second_coefficients;
                constraint.traction_first_coefficients = finite_region->traction_first_coefficients;
                constraint.traction_second_coefficients = finite_region->traction_second_coefficients;
                constraint.secondary_tangent_first_coefficients = finite_region->secondary_tangent_first_coefficients;
                constraint.secondary_tangent_second_coefficients = finite_region->secondary_tangent_second_coefficients;
                std::vector<double> primary_face_weight(_primary_contact_faces.at(constraint.contact).size(), 0.0);
                double friction_area = 0.0;
                for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
                    const SecondaryContactFace& secondary =
                        _secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
                    const Quad4FaceCoordinates secondary_current = current_face(secondary.nodes, secondary.coordinates);
                    const std::array<double, 2>& location =
                        abaqus_quad4_constraint_locations().at(sample.secondary_local_point);
                    const double local_area =
                        make_quad4_face_quadrature_point(secondary_current, location[0], location[1], 1.0)
                            .weighted_measure;
                    friction_area += local_area;
                    if (!sample.primary_faces.empty())
                        for (const std::size_t primary : sample.primary_faces)
                            primary_face_weight[primary] +=
                                local_area / static_cast<double>(sample.primary_faces.size());
                }
                if (!std::isfinite(friction_area) || !(friction_area > 0.0))
                    throw std::domain_error("HEX8 averaged friction constraint has a nonpositive current area");
                const double coefficient_scale = finite_region->area / friction_area;
                for (double& coefficient : constraint.gap_coefficients) coefficient *= coefficient_scale;
                for (double& coefficient : constraint.secondary_coefficients) coefficient *= coefficient_scale;
                for (std::size_t node = 0; node < constraint.nodes.size(); ++node)
                    for (std::size_t component = 0; component < 3; ++component) {
                        constraint.tangent_first_coefficients[node][component] *= coefficient_scale;
                        constraint.tangent_second_coefficients[node][component] *= coefficient_scale;
                        constraint.traction_first_coefficients[node][component] *= coefficient_scale;
                        constraint.traction_second_coefficients[node][component] *= coefficient_scale;
                    }
                for (std::size_t output = 0; output < constraint.secondary_output_nodes.size(); ++output)
                    for (std::size_t component = 0; component < 3; ++component) {
                        constraint.secondary_tangent_first_coefficients[output][component] *= coefficient_scale;
                        constraint.secondary_tangent_second_coefficients[output][component] *= coefficient_scale;
                    }
                constraint.area = friction_area;
                if (!primary_face_weight.empty())
                    constraint.primary_face = static_cast<std::size_t>(
                        std::max_element(primary_face_weight.begin(), primary_face_weight.end()) -
                        primary_face_weight.begin());
                CartesianPoint3 separation{};
                for (std::size_t node = 0; node < constraint.nodes.size(); ++node) {
                    separation.x += constraint.gap_coefficients[node] * constraint.reference_coordinates[node].x;
                    separation.y += constraint.gap_coefficients[node] * constraint.reference_coordinates[node].y;
                    separation.z += constraint.gap_coefficients[node] * constraint.reference_coordinates[node].z;
                }
                constraint.reference_gap = dot(separation, constraint.normal);
                continue;
            }
        }
        std::fill(constraint.gap_coefficients.begin(), constraint.gap_coefficients.end(), 0.0);
        std::fill(constraint.secondary_coefficients.begin(), constraint.secondary_coefficients.end(), 0.0);
        std::fill(constraint.tangent_first_coefficients.begin(), constraint.tangent_first_coefficients.end(),
            std::array<double, 3>{});
        std::fill(constraint.tangent_second_coefficients.begin(), constraint.tangent_second_coefficients.end(),
            std::array<double, 3>{});
        std::fill(constraint.traction_first_coefficients.begin(), constraint.traction_first_coefficients.end(),
            std::array<double, 3>{});
        std::fill(constraint.traction_second_coefficients.begin(), constraint.traction_second_coefficients.end(),
            std::array<double, 3>{});
        std::fill(constraint.secondary_tangent_first_coefficients.begin(),
            constraint.secondary_tangent_first_coefficients.end(), std::array<double, 3>{});
        std::fill(constraint.secondary_tangent_second_coefficients.begin(),
            constraint.secondary_tangent_second_coefficients.end(), std::array<double, 3>{});
        constraint.normal = {};
        constraint.tangent_first = {};
        constraint.area = 0.0;
        constraint.projected = true;
        const auto& primary_faces = _primary_contact_faces.at(constraint.contact);
        const auto& primary_current_faces = primary_current_faces_by_contact.at(constraint.contact);
        std::vector<double> primary_face_weight(primary_faces.size(), 0.0);
        for (AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
            const std::vector<std::size_t> previous_primary_faces = sample.primary_faces;
            const std::vector<AbaqusAveragedConstraint::FiniteSlidingSample::NormalPoint> previous_normal_points =
                sample.normal_points;
            sample.primary_faces.clear();
            sample.normal_points.clear();
            const SecondaryContactFace& secondary =
                _secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
            const Quad4FaceCoordinates secondary_current = current_face(secondary.nodes, secondary.coordinates);
            const std::array<double, 2>& location =
                abaqus_quad4_constraint_locations().at(sample.secondary_local_point);
            const Quad4FaceQuadraturePoint point =
                make_quad4_face_quadrature_point(secondary_current, location[0], location[1], 1.0);
            const double local_area = point.weighted_measure;
            constraint.area += local_area;
            for (std::size_t local_node = 0; local_node < 4; ++local_node) {
                const double coefficient = local_area * averaging[sample.secondary_local_point * 4 + local_node];
                const auto node =
                    std::find(constraint.nodes.begin(), constraint.nodes.end(), secondary.nodes[local_node]);
                if (node == constraint.nodes.end())
                    throw std::logic_error("HEX8 finite-sliding averaged constraint lost a secondary node");
                constraint.gap_coefficients[static_cast<std::size_t>(node - constraint.nodes.begin())] -= coefficient;
                const auto output = std::find(_secondary_boundaries.at(constraint.contact).boundary.nodes.begin(),
                    _secondary_boundaries.at(constraint.contact).boundary.nodes.end(),
                    _secondary_boundaries.at(constraint.contact)
                        .boundary.faces.at(sample.secondary_face)
                        .nodes[local_node]);
                if (output == _secondary_boundaries.at(constraint.contact).boundary.nodes.end())
                    throw std::logic_error("HEX8 finite-sliding averaged constraint lost an output node");
                const std::size_t output_index = static_cast<std::size_t>(
                    output - _secondary_boundaries.at(constraint.contact).boundary.nodes.begin());
                const auto stored = std::find(
                    constraint.secondary_output_nodes.begin(), constraint.secondary_output_nodes.end(), output_index);
                if (stored == constraint.secondary_output_nodes.end())
                    throw std::logic_error("HEX8 finite-sliding averaged constraint output support changed");
                constraint.secondary_coefficients[static_cast<std::size_t>(
                    stored - constraint.secondary_output_nodes.begin())] += coefficient;
            }
            const Quad4FaceQuadraturePoint normal_point = make_quad4_face_quadrature_point(
                secondary_current, (4.0 / 3.0) * location[0], (4.0 / 3.0) * location[1], 1.0);
            const Quad4FaceQuadraturePoint reference_point =
                make_quad4_face_quadrature_point(secondary.coordinates, location[0], location[1], 1.0);
            const CartesianPoint3 area_vector = cross(normal_point.tangent_xi, normal_point.tangent_eta);
            const double normal_measure = std::sqrt(dot(area_vector, area_vector));
            if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
                throw std::domain_error("HEX8 finite-sliding averaged constraint has an undefined current normal");
            constraint.normal.x += local_area * sample.normal_orientation * area_vector.x / normal_measure;
            constraint.normal.y += local_area * sample.normal_orientation * area_vector.y / normal_measure;
            constraint.normal.z += local_area * sample.normal_orientation * area_vector.z / normal_measure;
            const CartesianPoint3 local_normal = {sample.normal_orientation * area_vector.x / normal_measure,
                sample.normal_orientation * area_vector.y / normal_measure,
                sample.normal_orientation * area_vector.z / normal_measure};
            CartesianPoint3 tangent = point.tangent_xi;
            const double tangent_measure = std::sqrt(dot(tangent, tangent));
            if (!std::isfinite(tangent_measure) || !(tangent_measure > 0.0))
                throw std::domain_error("HEX8 finite-sliding averaged constraint has an undefined current tangent");
            const CartesianPoint3 local_first = {sample.tangent_orientation * tangent.x / tangent_measure,
                                      sample.tangent_orientation * tangent.y / tangent_measure,
                                      sample.tangent_orientation * tangent.z / tangent_measure},
                                  local_second = cross(local_normal, local_first);
            const SurfacePullbackCovectors pullback = surface_pullback_covectors(point, reference_point,
                sample.normal_orientation, sample.tangent_orientation, local_normal, local_first);
            constraint.tangent_first.x += local_area * local_first.x;
            constraint.tangent_first.y += local_area * local_first.y;
            constraint.tangent_first.z += local_area * local_first.z;
            for (std::size_t local_node = 0; local_node < 4; ++local_node) {
                const double coefficient = local_area * averaging[sample.secondary_local_point * 4 + local_node];
                const auto node =
                    std::find(constraint.nodes.begin(), constraint.nodes.end(), secondary.nodes[local_node]);
                const std::size_t node_index = static_cast<std::size_t>(node - constraint.nodes.begin());
                const auto output = std::find(_secondary_boundaries.at(constraint.contact).boundary.nodes.begin(),
                    _secondary_boundaries.at(constraint.contact).boundary.nodes.end(),
                    _secondary_boundaries.at(constraint.contact)
                        .boundary.faces.at(sample.secondary_face)
                        .nodes[local_node]);
                const std::size_t output_index = static_cast<std::size_t>(
                    output - _secondary_boundaries.at(constraint.contact).boundary.nodes.begin());
                const auto stored = std::find(
                    constraint.secondary_output_nodes.begin(), constraint.secondary_output_nodes.end(), output_index);
                const std::size_t stored_index =
                    static_cast<std::size_t>(stored - constraint.secondary_output_nodes.begin());
                for (std::size_t component = 0; component < 3; ++component) {
                    const double pullback_first = component == 0   ? pullback.first.x
                                                  : component == 1 ? pullback.first.y
                                                                   : pullback.first.z,
                                 pullback_second = component == 0   ? pullback.second.x
                                                   : component == 1 ? pullback.second.y
                                                                    : pullback.second.z,
                                 first = component == 0   ? local_first.x
                                         : component == 1 ? local_first.y
                                                          : local_first.z,
                                 second = component == 0   ? local_second.x
                                          : component == 1 ? local_second.y
                                                           : local_second.z;
                    constraint.tangent_first_coefficients[node_index][component] -= coefficient * pullback_first;
                    constraint.tangent_second_coefficients[node_index][component] -= coefficient * pullback_second;
                    constraint.traction_first_coefficients[node_index][component] -= coefficient * first;
                    constraint.traction_second_coefficients[node_index][component] -= coefficient * second;
                    constraint.secondary_tangent_first_coefficients[stored_index][component] += coefficient * first;
                    constraint.secondary_tangent_second_coefficients[stored_index][component] += coefficient * second;
                }
            }

            struct Projection final {
                std::size_t primary;
                Quad4ReferenceProjectionValue value;
                bool interior;
            };

            std::vector<Projection> projections;
            double minimum_projection_distance = std::numeric_limits<double>::infinity();
            std::vector<unsigned char> considered(primary_faces.size(), 0U);
            const auto consider_projection = [&](std::size_t primary) {
                if (considered[primary] != 0U) return;
                considered[primary] = 1U;
                const Quad4ReferenceProjectionValue value = compute_quad4_reference_projection(
                    secondary_current, primary_current_faces[primary], point.shape, 1.0);
                if (!value.projected) return;
                const double minimum_shape = *std::min_element(value.primary_shape.begin(), value.primary_shape.end());
                projections.push_back({primary, value, minimum_shape > boundary_tolerance});
                minimum_projection_distance = std::min(minimum_projection_distance, std::abs(value.gap));
            };
            for (const std::size_t primary : previous_primary_faces)
                if (primary < primary_faces.size()) consider_projection(primary);
            if (primary_faces.size() > spatial_detail::contact_search_tree_minimum_items) {
                std::array<double, 3> search_point{};
                for (std::size_t node = 0; node < 4; ++node) {
                    search_point[0] += point.shape[node] * secondary_current[node].x;
                    search_point[1] += point.shape[node] * secondary_current[node].y;
                    search_point[2] += point.shape[node] * secondary_current[node].z;
                }
                _contact_search_trees[constraint.contact].begin_query(search_point, _contact_search_query);
                std::size_t primary = 0;
                double candidate_distance = minimum_projection_distance;
                while (_contact_search_trees[constraint.contact].next_candidate(
                    _contact_search_query, candidate_distance, primary)) {
                    consider_projection(primary);
                    candidate_distance = minimum_projection_distance;
                }
            } else
                for (std::size_t primary = 0; primary < primary_faces.size(); ++primary) consider_projection(primary);
            if (projections.empty()) {
                constraint.projected = false;
                continue;
            }
            const bool has_interior = std::any_of(
                projections.begin(), projections.end(), [](const Projection& value) { return value.interior; });
            if (has_interior && projections.size() > 1) {
                const auto closest = std::min_element(
                    projections.begin(), projections.end(), [](const Projection& a, const Projection& b) {
                        return std::abs(a.value.gap) < std::abs(b.value.gap);
                    });
                const Projection selected = *closest;
                projections.assign(1, selected);
            }
            const double projection_weight = local_area / static_cast<double>(projections.size());
            for (const Projection& projection : projections) {
                sample.primary_faces.push_back(projection.primary);
                const PrimaryContactFace& face = primary_faces[projection.primary];
                primary_face_weight[projection.primary] += projection_weight;
                for (std::size_t local_node = 0; local_node < 4; ++local_node) {
                    const auto node =
                        std::find(constraint.nodes.begin(), constraint.nodes.end(), face.nodes[local_node]);
                    if (node == constraint.nodes.end())
                        throw std::logic_error("HEX8 finite-sliding averaged constraint lost a primary node");
                    constraint.gap_coefficients[static_cast<std::size_t>(node - constraint.nodes.begin())] +=
                        projection_weight * projection.value.primary_shape[local_node];
                    const std::size_t node_index = static_cast<std::size_t>(node - constraint.nodes.begin());
                    for (std::size_t component = 0; component < 3; ++component) {
                        const double coefficient = projection_weight * projection.value.primary_shape[local_node],
                                     pullback_first = component == 0   ? pullback.first.x
                                                      : component == 1 ? pullback.first.y
                                                                       : pullback.first.z,
                                     pullback_second = component == 0   ? pullback.second.x
                                                       : component == 1 ? pullback.second.y
                                                                        : pullback.second.z,
                                     first = component == 0   ? local_first.x
                                             : component == 1 ? local_first.y
                                                              : local_first.z,
                                     second = component == 0   ? local_second.x
                                              : component == 1 ? local_second.y
                                                               : local_second.z;
                        constraint.tangent_first_coefficients[node_index][component] += coefficient * pullback_first;
                        constraint.tangent_second_coefficients[node_index][component] += coefficient * pullback_second;
                        constraint.traction_first_coefficients[node_index][component] += coefficient * first;
                        constraint.traction_second_coefficients[node_index][component] += coefficient * second;
                    }
                }
            }
            if (constraint.finite_region_normal) {
                constexpr std::array<double, 2> integration_points = {-0.5773502691896258, 0.5773502691896258};
                constexpr std::array<double, 2> integration_weights = {1.0, 1.0};
                const double xi_center = location[0], eta_center = location[1];
                const auto shared_edge = [&](std::size_t primary, std::size_t first, std::size_t second) {
                    constexpr std::array<std::array<std::size_t, 2>, 4> edges = {
                        {{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}}};
                    for (std::size_t edge = 0; edge < edges.size(); ++edge)
                        if ((edges[edge][0] == first && edges[edge][1] == second) ||
                            (edges[edge][0] == second && edges[edge][1] == first))
                            return primary_faces[primary].shared_edges[edge];
                    throw std::logic_error("HEX8 finite-region normal search requested an invalid primary edge");
                };
                const auto owner = [&](double xi, double eta, std::size_t cached_primary) {
                    const Quad4FaceQuadraturePoint integration_point =
                        make_quad4_face_quadrature_point(secondary_current, xi, eta, 1.0);
                    double minimum_distance = std::numeric_limits<double>::infinity();
                    std::size_t selected = std::numeric_limits<std::size_t>::max();
                    const auto consider = [&](std::size_t primary) {
                        const Quad4ReferenceProjectionValue unbounded = compute_quad4_reference_projection(
                            secondary_current, primary_current_faces[primary], integration_point.shape, 1.0);
                        Quad4ReferenceProjectionValue value = unbounded;
                        if (!value.projected) {
                            const bool xi_supported = (unbounded.xi >= -1.0 && unbounded.xi <= 1.0) ||
                                                      (unbounded.xi < -1.0 && shared_edge(primary, 0, 3)) ||
                                                      (unbounded.xi > 1.0 && shared_edge(primary, 1, 2)),
                                       eta_supported = (unbounded.eta >= -1.0 && unbounded.eta <= 1.0) ||
                                                       (unbounded.eta < -1.0 && shared_edge(primary, 0, 1)) ||
                                                       (unbounded.eta > 1.0 && shared_edge(primary, 3, 2));
                            if (!xi_supported || !eta_supported) return;
                            value = compute_quad4_reference_closest_projection(
                                secondary_current, primary_current_faces[primary], integration_point.shape, 1.0);
                        }
                        if (!value.projected) return;
                        const double distance = value.distance;
                        if (distance < minimum_distance || (distance == minimum_distance && primary < selected)) {
                            minimum_distance = distance;
                            selected = primary;
                        }
                    };
                    if (cached_primary < primary_faces.size()) consider(cached_primary);
                    if (primary_faces.size() > spatial_detail::contact_search_tree_minimum_items) {
                        std::array<double, 3> search_point{};
                        for (std::size_t node = 0; node < 4; ++node) {
                            search_point[0] += integration_point.shape[node] * secondary_current[node].x;
                            search_point[1] += integration_point.shape[node] * secondary_current[node].y;
                            search_point[2] += integration_point.shape[node] * secondary_current[node].z;
                        }
                        _contact_search_trees[constraint.contact].begin_query(search_point, _contact_search_query);
                        std::size_t primary = 0;
                        double candidate_distance = minimum_distance;
                        while (_contact_search_trees[constraint.contact].next_candidate(
                            _contact_search_query, candidate_distance, primary)) {
                            if (primary != cached_primary) consider(primary);
                            candidate_distance = minimum_distance;
                        }
                    } else
                        for (std::size_t primary = 0; primary < primary_faces.size(); ++primary)
                            if (primary != cached_primary) consider(primary);
                    return selected;
                };
                constexpr double region_half_width = 0.5;
                std::size_t normal_point_index = 0;
                for (std::size_t eta_point = 0; eta_point < integration_points.size(); ++eta_point)
                    for (std::size_t xi_point = 0; xi_point < integration_points.size(); ++xi_point) {
                        const double xi = xi_center + region_half_width * integration_points[xi_point],
                                     eta = eta_center + region_half_width * integration_points[eta_point],
                                     weight = region_half_width * region_half_width * integration_weights[xi_point] *
                                              integration_weights[eta_point];
                        const std::size_t cached_primary = normal_point_index < previous_normal_points.size()
                                                               ? previous_normal_points[normal_point_index].primary_face
                                                               : std::numeric_limits<std::size_t>::max();
                        const std::size_t primary = owner(xi, eta, cached_primary);
                        ++normal_point_index;
                        if (primary == std::numeric_limits<std::size_t>::max()) {
                            constraint.projected = false;
                            continue;
                        }
                        sample.normal_points.push_back({xi, eta, weight, primary});
                        primary_face_weight[primary] += weight;
                    }
            }
        }
        if (constraint.finite_region_normal) {
            constraint.area = 0.0;
            for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
                const SecondaryContactFace& secondary =
                    _secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
                const Quad4FaceCoordinates secondary_current = current_face(secondary.nodes, secondary.coordinates);
                for (const AbaqusAveragedConstraint::FiniteSlidingSample::NormalPoint& normal_point :
                    sample.normal_points)
                    constraint.area += make_quad4_face_quadrature_point(
                        secondary_current, normal_point.xi, normal_point.eta, normal_point.weight)
                                           .weighted_measure;
            }
        }
        if (!constraint.projected) continue;
        if (!std::isfinite(constraint.area) || !(constraint.area > 0.0))
            throw std::domain_error("HEX8 finite-sliding averaged constraint has a nonpositive current area");
        const double normal_measure = std::sqrt(dot(constraint.normal, constraint.normal));
        if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
            throw std::domain_error("HEX8 finite-sliding averaged constraint has an undefined averaged current normal");
        constraint.normal.x /= normal_measure;
        constraint.normal.y /= normal_measure;
        constraint.normal.z /= normal_measure;
        const double tangent_normal_component = dot(constraint.tangent_first, constraint.normal);
        constraint.tangent_first.x -= tangent_normal_component * constraint.normal.x;
        constraint.tangent_first.y -= tangent_normal_component * constraint.normal.y;
        constraint.tangent_first.z -= tangent_normal_component * constraint.normal.z;
        const double tangent_measure = std::sqrt(dot(constraint.tangent_first, constraint.tangent_first));
        if (!std::isfinite(tangent_measure) || !(tangent_measure > 0.0))
            throw std::domain_error(
                "HEX8 finite-sliding averaged constraint has an undefined averaged current tangent");
        constraint.tangent_first.x /= tangent_measure;
        constraint.tangent_first.y /= tangent_measure;
        constraint.tangent_first.z /= tangent_measure;
        for (double& coefficient : constraint.gap_coefficients) coefficient /= constraint.area;
        for (double& coefficient : constraint.secondary_coefficients) coefficient /= constraint.area;
        for (std::size_t node = 0; node < constraint.nodes.size(); ++node)
            for (std::size_t component = 0; component < 3; ++component) {
                constraint.tangent_first_coefficients[node][component] /= constraint.area;
                constraint.tangent_second_coefficients[node][component] /= constraint.area;
                constraint.traction_first_coefficients[node][component] /= constraint.area;
                constraint.traction_second_coefficients[node][component] /= constraint.area;
            }
        for (std::size_t output = 0; output < constraint.secondary_output_nodes.size(); ++output)
            for (std::size_t component = 0; component < 3; ++component) {
                constraint.secondary_tangent_first_coefficients[output][component] /= constraint.area;
                constraint.secondary_tangent_second_coefficients[output][component] /= constraint.area;
            }
        constraint.primary_face = static_cast<std::size_t>(
            std::max_element(primary_face_weight.begin(), primary_face_weight.end()) - primary_face_weight.begin());
        CartesianPoint3 separation{};
        for (std::size_t node = 0; node < constraint.nodes.size(); ++node) {
            separation.x += constraint.gap_coefficients[node] * constraint.reference_coordinates[node].x;
            separation.y += constraint.gap_coefficients[node] * constraint.reference_coordinates[node].y;
            separation.z += constraint.gap_coefficients[node] * constraint.reference_coordinates[node].z;
        }
        constraint.reference_gap = dot(separation, constraint.normal);
    }
    const auto set_active_support = [](AbaqusAveragedConstraint& constraint, std::vector<std::size_t> required_nodes) {
        std::sort(required_nodes.begin(), required_nodes.end());
        required_nodes.erase(std::unique(required_nodes.begin(), required_nodes.end()), required_nodes.end());
        constraint.active_nodes.clear();
        constraint.active_node_indices.clear();
        for (std::size_t stored = 0; stored < constraint.nodes.size(); ++stored)
            if (std::binary_search(required_nodes.begin(), required_nodes.end(), constraint.nodes[stored])) {
                constraint.active_nodes.push_back(constraint.nodes[stored]);
                constraint.active_node_indices.push_back(stored);
            }
        if (constraint.active_nodes.empty())
            throw std::logic_error("HEX8 finite-sliding averaged constraint has no active nodes");
    };
    for (AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        if (!constraint.finite_sliding) continue;
        std::vector<std::size_t> required_nodes;
        for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
            const SecondaryContactFace& secondary =
                _secondary_contact_faces.at(constraint.contact).at(sample.secondary_face);
            required_nodes.insert(required_nodes.end(), secondary.nodes.begin(), secondary.nodes.end());
            for (const std::size_t primary_face : sample.primary_faces) {
                const PrimaryContactFace& primary = _primary_contact_faces.at(constraint.contact).at(primary_face);
                required_nodes.insert(required_nodes.end(), primary.nodes.begin(), primary.nodes.end());
            }
            for (const AbaqusAveragedConstraint::FiniteSlidingSample::NormalPoint& normal_point :
                sample.normal_points) {
                const PrimaryContactFace& primary =
                    _primary_contact_faces.at(constraint.contact).at(normal_point.primary_face);
                required_nodes.insert(required_nodes.end(), primary.nodes.begin(), primary.nodes.end());
            }
        }
        set_active_support(constraint, std::move(required_nodes));
    }
    for (AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        if (!constraint.friction_only) continue;
        const auto finite_region = std::find_if(_abaqus_averaged_constraints.begin(),
            _abaqus_averaged_constraints.end(), [&constraint](const AbaqusAveragedConstraint& value) {
                return value.contact == constraint.contact && value.secondary == constraint.secondary &&
                       value.finite_region_normal;
            });
        if (finite_region == _abaqus_averaged_constraints.end()) continue;
        std::vector<std::size_t> required_nodes = constraint.active_nodes;
        required_nodes.insert(
            required_nodes.end(), finite_region->active_nodes.begin(), finite_region->active_nodes.end());
        set_active_support(constraint, std::move(required_nodes));
    }
}

bool SpatialAssembly::summarize_averaged_contact(std::size_t contact_value, const std::vector<double>& state,
    std::vector<CartesianContactNodeSummary>& summaries) const {
    refresh_finite_averaged_constraints(state);
    const bool is_committed_state = state == _committed_contact_solution;
    const bool averaged = std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
        [contact_value](const AbaqusAveragedConstraint& value) { return value.contact == contact_value; });
    if (!averaged) return false;
    const bool separately_averaged_friction = std::any_of(_abaqus_averaged_constraints.begin(),
        _abaqus_averaged_constraints.end(), [contact_value](const AbaqusAveragedConstraint& value) {
            return value.contact == contact_value && value.friction_only;
        });
    std::vector<std::size_t> dofs;
    std::vector<std::array<double, 3>> finite_region_normal_area(summaries.size());
    bool recover_finite_sliding_nodal_tractions = false;
    bool finite_region_normals = false;
    for (const AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        if (constraint.contact != contact_value) continue;
        recover_finite_sliding_nodal_tractions = recover_finite_sliding_nodal_tractions || constraint.finite_sliding;
        averaged_constraint_dofs(constraint, dofs);
        std::vector<double> local_state(dofs.size()), committed_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local) {
            local_state[local] = state.at(dofs[local]);
            committed_state[local] = _committed_contact_solution.at(dofs[local]);
        }
        std::vector<double> finite_region_residual;
        const AbaqusAveragedConstraintValue value =
            constraint.finite_region_normal
                ? finite_region_normal_value(constraint, constraint.active_nodes, local_state, &finite_region_residual)
                : averaged_constraint_value(constraint, local_state, committed_state,
                      _contact_histories[constraint.contact][constraint.history]);
        CartesianContactNodeSummary& predominant = summaries.at(constraint.secondary);
        predominant.projected = constraint.projected;
        if (!constraint.friction_only) {
            predominant.primary_face = constraint.primary_face;
            predominant.gap = value.gap;
            predominant.constraint_pressure = value.pressure;
            predominant.pressure = value.pressure;
            predominant.tributary_area = constraint.area;
        }
        if (constraint.finite_region_normal) {
            finite_region_normals = true;
            const ResolvedBoundary& secondary = _secondary_boundaries.at(constraint.contact);
            for (std::size_t entry = 0; entry < constraint.secondary_output_nodes.size(); ++entry) {
                const std::size_t output_index = constraint.secondary_output_nodes[entry];
                const std::size_t global = global_node(secondary.region, secondary.boundary.nodes.at(output_index));
                const auto stored = std::find(constraint.active_nodes.begin(), constraint.active_nodes.end(), global);
                if (stored == constraint.active_nodes.end())
                    throw std::logic_error("HEX8 finite-region normal output support changed");
                const std::size_t stored_index = static_cast<std::size_t>(stored - constraint.active_nodes.begin());
                CartesianContactNodeSummary& output = summaries.at(output_index);
                output.projected = output.projected || constraint.projected;
                output.primary_face = constraint.primary_face;
                output.contact_force += constraint.secondary_coefficients[entry] * value.force;
                for (std::size_t component = 0; component < 3; ++component)
                    output.normal_contact_force[component] +=
                        finite_region_residual[component * constraint.active_nodes.size() + stored_index];
                if (value.pressure > 0.0)
                    for (std::size_t component = 0; component < 3; ++component)
                        finite_region_normal_area[output_index][component] +=
                            finite_region_residual[component * constraint.active_nodes.size() + stored_index] /
                            value.pressure;
            }
            continue;
        }
        if (!separately_averaged_friction || constraint.friction_only) {
            predominant.tangential_slip = value.tangential_slip;
            predominant.elastic_tangential_slip = value.elastic_tangential_slip;
            predominant.sliding =
                is_committed_state ? _contact_histories[constraint.contact][constraint.history].sliding : value.sliding;
        }
        for (std::size_t entry = 0; entry < constraint.secondary_output_nodes.size(); ++entry) {
            CartesianContactNodeSummary& output = summaries.at(constraint.secondary_output_nodes[entry]);
            output.projected = output.projected || constraint.projected;
            output.primary_face = constraint.primary_face;
            if (!constraint.friction_only)
                output.contact_force += constraint.secondary_coefficients[entry] * value.force;
            output.tangential_force += constraint.secondary_coefficients[entry] * value.tangential_force;
            for (std::size_t component = 0; component < 3; ++component) {
                const double normal = component == 0   ? constraint.normal.x
                                      : component == 1 ? constraint.normal.y
                                                       : constraint.normal.z;
                if (!constraint.friction_only)
                    output.normal_contact_force[component] +=
                        constraint.secondary_normal_coefficients.empty()
                            ? constraint.secondary_coefficients[entry] * value.force * normal
                            : constraint.secondary_normal_coefficients[entry][component] * value.force;
                output.tangential_contact_force[component] +=
                    constraint.friction_only
                        ? constraint.area * (constraint.secondary_tangent_first_coefficients[entry][component] *
                                                    value.tangential_traction_components[0] +
                                                constraint.secondary_tangent_second_coefficients[entry][component] *
                                                    value.tangential_traction_components[1])
                        : constraint.secondary_coefficients[entry] * constraint.area *
                              value.tangential_traction[component];
            }
        }
    }
    if (recover_finite_sliding_nodal_tractions && !finite_region_normals &&
        std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
            [contact_value](const AbaqusAveragedConstraint& constraint) {
                return constraint.contact == contact_value && !constraint.friction_only;
            }))
        for (CartesianContactNodeSummary& summary : summaries)
            if (summary.tributary_area > 0.0) {
                if (!_uses_hex20) summary.pressure = summary.contact_force / summary.tributary_area;
                summary.tangential_traction = summary.tangential_force / summary.tributary_area;
            }
    if (_uses_hex20 && recover_finite_sliding_nodal_tractions && !finite_region_normals) {
        const ResolvedHex20Boundary& secondary = _hex20_secondary_boundaries.at(contact_value);
        const Matrix8& recovery = abaqus_quad8_pressure_recovery();
        std::vector<double> recovered_pressure(summaries.size());
        std::vector<std::size_t> recovery_count(summaries.size());
        std::vector<std::size_t> component_parent(summaries.size());
        std::iota(component_parent.begin(), component_parent.end(), 0);
        const auto component_root = [&component_parent](std::size_t node) {
            while (component_parent[node] != node) {
                component_parent[node] = component_parent[component_parent[node]];
                node = component_parent[node];
            }
            return node;
        };
        for (const Quad8FaceElement& face : secondary.boundary.faces) {
            std::array<std::size_t, 8> output_nodes{};
            for (std::size_t local_node = 0; local_node < face.nodes.size(); ++local_node) {
                const auto found = std::find(secondary.boundary.displacement_nodes.begin(),
                    secondary.boundary.displacement_nodes.end(), face.nodes[local_node]);
                if (found == secondary.boundary.displacement_nodes.end())
                    throw std::logic_error("HEX20 contact-pressure recovery node mapping failed");
                const std::size_t output_node =
                    static_cast<std::size_t>(found - secondary.boundary.displacement_nodes.begin());
                output_nodes[local_node] = output_node;
            }
            for (std::size_t local_node = 1; local_node < output_nodes.size(); ++local_node) {
                const std::size_t first_root = component_root(output_nodes[0]);
                const std::size_t next_root = component_root(output_nodes[local_node]);
                if (first_root != next_root) component_parent[next_root] = first_root;
            }
            for (std::size_t local_node = 0; local_node < output_nodes.size(); ++local_node) {
                const std::size_t output_node = output_nodes[local_node];
                for (std::size_t input_node = 0; input_node < output_nodes.size(); ++input_node) {
                    const double pressure = summaries[output_nodes[input_node]].constraint_pressure;
                    recovered_pressure[output_node] += recovery[local_node * 8 + input_node] * pressure;
                }
                ++recovery_count[output_node];
            }
        }
        std::vector<double> component_minimum(summaries.size(), std::numeric_limits<double>::infinity());
        std::vector<double> component_maximum(summaries.size(), -std::numeric_limits<double>::infinity());
        for (std::size_t node = 0; node < summaries.size(); ++node)
            if (recovery_count[node] != 0) {
                const std::size_t root = component_root(node);
                component_minimum[root] = std::min(component_minimum[root], summaries[node].constraint_pressure);
                component_maximum[root] = std::max(component_maximum[root], summaries[node].constraint_pressure);
            }
        for (std::size_t node = 0; node < summaries.size(); ++node)
            if (recovery_count[node] != 0) {
                // H20.42 shows that face values are averaged arithmetically at shared nodes. Component-wide bounds
                // reproduce partial-contact probes and provide a conservative surrogate for Abaqus's unidentified
                // nonlinear branch; all 128 H20.42 states remain below one percent in each accepted relative metric.
                const std::size_t root = component_root(node);
                const double average = recovered_pressure[node] / static_cast<double>(recovery_count[node]);
                summaries[node].pressure = std::clamp(average, component_minimum[root], component_maximum[root]);
            }
    }
    if (finite_region_normals)
        for (std::size_t node = 0; node < summaries.size(); ++node) {
            CartesianContactNodeSummary& summary = summaries[node];
            const double normal_force = std::hypot(
                summary.normal_contact_force[0], summary.normal_contact_force[1], summary.normal_contact_force[2]);
            summary.tributary_area = std::hypot(finite_region_normal_area[node][0], finite_region_normal_area[node][1],
                finite_region_normal_area[node][2]);
            if (summary.tributary_area > 0.0) {
                summary.pressure = normal_force / summary.tributary_area;
                summary.constraint_pressure = summary.pressure;
                summary.tangential_traction = summary.tangential_force / summary.tributary_area;
            }
        }
    return std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
        [contact_value](const AbaqusAveragedConstraint& constraint) {
            return constraint.contact == contact_value && !constraint.friction_only;
        });
}

SpatialAssembly::FiniteRegionPartitionSummary SpatialAssembly::finite_region_partition_summary(
    std::size_t contact_index, const std::vector<double>& state) const {
    refresh_finite_averaged_constraints(state);
    FiniteRegionPartitionSummary result{};
    std::set<std::size_t> active_primary_faces;
    for (const AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        if (constraint.contact != contact_index || !constraint.finite_region_normal) continue;
        ++result.constraint_count;
        result.all_projected = result.all_projected && constraint.projected;
        std::set<std::size_t> constraint_primary_faces;
        for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples)
            for (const AbaqusAveragedConstraint::FiniteSlidingSample::NormalPoint& point : sample.normal_points) {
                ++result.integration_point_count;
                result.maximum_owners_per_integration_point = 1;
                active_primary_faces.insert(point.primary_face);
                constraint_primary_faces.insert(point.primary_face);
            }
        if (constraint_primary_faces.size() > 1) ++result.cross_face_constraint_count;
    }
    result.active_primary_face_count = active_primary_faces.size();
    return result;
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
        if (definition.mechanical_discretization == MechanicalContactDiscretization::automatic)
            definition.mechanical_discretization = MechanicalContactDiscretization::node_to_surface;
        const bool surface_to_surface =
            definition.mechanical_discretization == MechanicalContactDiscretization::surface_to_surface;
        const bool finite_sliding =
            surface_to_surface && definition.mechanical_sliding == MechanicalContactSliding::finite;
        bool finite_averaged = definition.mechanical && finite_sliding && definition.friction_coefficient == 0.0;
        if (definition.friction_slip_tolerance > 0.0 && !surface_to_surface)
            throw std::invalid_argument("slip_tolerance requires surface_to_surface contact: " + definition.name);
        if (definition.quad8_nodal_area_rule != Quad8NodalAreaRule::positive_lumped)
            throw std::invalid_argument(
                "quad8_nodal_area_rule = consistent_shape is supported only for HEX20 contact comparisons: " +
                definition.name);
        ResolvedBoundary primary = resolve_boundary(source_mesh, definition.primary),
                         secondary = resolve_boundary(source_mesh, definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("Three-dimensional self-contact is not supported: " + definition.name);
        if (surface_to_surface && !finite_sliding &&
            (_definition.regions[primary.region].strain_formulation != StrainFormulation::small ||
                _definition.regions[secondary.region].strain_formulation != StrainFormulation::small))
            throw std::invalid_argument(
                "HEX8 small-sliding surface_to_surface contact currently requires small-strain regions: " +
                definition.name);
        const Hex8RegionMesh& primary_mesh = _meshes[primary.region];
        const Hex8RegionMesh& secondary_mesh = _meshes[secondary.region];
        bool noncoplanar_secondary = false;
        if (secondary.boundary.faces.size() > 1) {
            const Quad4FaceQuadraturePoint reference_normal_point = make_quad4_face_quadrature_point(
                face_coordinates(secondary_mesh, secondary.boundary.faces.front()), 0.0, 0.0, 1.0);
            const CartesianPoint3 reference_normal =
                cross(reference_normal_point.tangent_xi, reference_normal_point.tangent_eta);
            const double reference_measure = std::sqrt(dot(reference_normal, reference_normal));
            for (std::size_t face = 1; face < secondary.boundary.faces.size(); ++face) {
                const Quad4FaceQuadraturePoint normal_point = make_quad4_face_quadrature_point(
                    face_coordinates(secondary_mesh, secondary.boundary.faces[face]), 0.0, 0.0, 1.0);
                const CartesianPoint3 normal = cross(normal_point.tangent_xi, normal_point.tangent_eta);
                const double measure = std::sqrt(dot(normal, normal));
                const double cosine = std::abs(dot(reference_normal, normal) / (reference_measure * measure));
                if (cosine < 1.0 - 1.0e-10) {
                    noncoplanar_secondary = true;
                    break;
                }
            }
        }
        const bool finite_region_normal = definition.mechanical && finite_sliding && noncoplanar_secondary;
        if (finite_region_normal) finite_averaged = true;
        const bool finite_averaged_friction = finite_region_normal && definition.friction_coefficient > 0.0;
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
        _thermal_properties.push_back(
            {definition.thermal ? definition.gap_conductivity : 1.0, definition.thermal ? definition.minimum_gap : 1.0,
                definition.thermal ? definition.gap_heat_conductance_law : GapHeatConductanceLaw::gas_gap,
                definition.thermal ? definition.gap_conductance : 0.0,
                definition.thermal ? definition.gap_conductance_clearance_derivative : 0.0,
                definition.thermal ? definition.gap_conductance_pressure_derivative : 0.0,
                definition.thermal ? definition.gap_conductance_temperature_derivative : 0.0,
                definition.thermal ? definition.gap_conductance_reference_temperature : 0.0,
                definition.mechanical ? definition.penalty : 0.0});
        const double slip_tolerance =
            definition.mechanical && definition.friction_coefficient > 0.0 && surface_to_surface
                ? (definition.friction_slip_tolerance > 0.0 ? definition.friction_slip_tolerance : 5.0e-3)
                : 0.0;
        double reference_secondary_area = 0.0;
        if (slip_tolerance > 0.0)
            for (const Quad4FaceElement& face : secondary.boundary.faces)
                for (const Quad4FaceQuadraturePoint& point :
                    make_quad4_face_geometry(face_coordinates(secondary_mesh, face)).points)
                    reference_secondary_area += point.weighted_measure;
        const double maximum_elastic_slip =
            slip_tolerance > 0.0 ? slip_tolerance * std::sqrt(reference_secondary_area /
                                                              static_cast<double>(secondary.boundary.faces.size()))
                                 : 0.0;
        if (!std::isfinite(maximum_elastic_slip) || maximum_elastic_slip < 0.0)
            throw std::invalid_argument(
                "Three-dimensional contact slip_tolerance gives an invalid elastic slip: " + definition.name);
        _mechanical_properties.push_back({definition.mechanical ? definition.penalty : 1.0,
            definition.mechanical ? definition.friction_coefficient : 0.0, false, maximum_elastic_slip});
        std::vector<PrimaryContactFace> primary_faces;
        primary_faces.reserve(primary.boundary.faces.size());
        for (const Quad4FaceElement& primary_face : primary.boundary.faces) {
            const Quad4FaceCoordinates coordinates = face_coordinates(primary_mesh, primary_face);
            std::array<std::size_t, 4> nodes{};
            for (std::size_t node = 0; node < nodes.size(); ++node)
                nodes[node] = global_node(primary.region, primary_face.nodes[node]);
            primary_faces.push_back({nodes, coordinates, element_centroid(primary_mesh, primary_face.parent_element)});
        }
        constexpr std::array<std::array<std::size_t, 2>, 4> quad4_edges = {{{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}}};
        for (std::size_t face = 0; face < primary_faces.size(); ++face)
            for (std::size_t edge = 0; edge < quad4_edges.size(); ++edge) {
                const std::size_t first = primary_faces[face].nodes[quad4_edges[edge][0]],
                                  second = primary_faces[face].nodes[quad4_edges[edge][1]];
                primary_faces[face].shared_edges[edge] =
                    std::any_of(primary_faces.begin(), primary_faces.end(), [&](const PrimaryContactFace& other) {
                        if (&other == &primary_faces[face]) return false;
                        return std::find(other.nodes.begin(), other.nodes.end(), first) != other.nodes.end() &&
                               std::find(other.nodes.begin(), other.nodes.end(), second) != other.nodes.end();
                    });
            }
        if (finite_averaged && !finite_region_normal && definition.friction_coefficient == 0.0 &&
            primary_faces.size() > 1) {
            const Quad4FaceQuadraturePoint reference_normal_point =
                make_quad4_face_quadrature_point(primary_faces.front().coordinates, 0.0, 0.0, 1.0);
            const CartesianPoint3 reference_normal =
                cross(reference_normal_point.tangent_xi, reference_normal_point.tangent_eta);
            const double reference_measure = std::sqrt(dot(reference_normal, reference_normal));
            for (std::size_t face = 1; face < primary_faces.size(); ++face) {
                const Quad4FaceQuadraturePoint normal_point =
                    make_quad4_face_quadrature_point(primary_faces[face].coordinates, 0.0, 0.0, 1.0);
                const CartesianPoint3 normal = cross(normal_point.tangent_xi, normal_point.tangent_eta);
                const double measure = std::sqrt(dot(normal, normal));
                const double cosine = std::abs(dot(reference_normal, normal) / (reference_measure * measure));
                if (cosine < 1.0 - 1.0e-10) {
                    finite_averaged = false;
                    break;
                }
            }
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
            std::array<Quad4FaceQuadraturePoint, 4> contact_points{};
            for (std::size_t point = 0; point < contact_points.size(); ++point) {
                const std::array<double, 2>& location = abaqus_quad4_constraint_locations()[point];
                contact_points[point] =
                    make_quad4_face_quadrature_point(secondary_coordinates, location[0], location[1], 1.0);
            }
            secondary_faces.push_back(
                {secondary_nodes, secondary_coordinates, secondary_geometry, secondary_parent, contact_points});
            if (definition.thermal) thermal_point_count += secondary_geometry.points.size();
            if (definition.mechanical && !surface_to_surface)
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
            if (definition.mechanical && finite_sliding && !finite_averaged)
                for (std::size_t point = 0; point < contact_points.size(); ++point) {
                    const std::size_t history = _contact_histories[contact_value].size();
                    _contact_histories[contact_value].emplace_back();
                    _mechanical_points.push_back({contact_value, history, secondary_face_index, point});
                }
        }
        if (!finite_sliding || finite_averaged)
            _contact_histories[contact_value].resize(secondary.boundary.nodes.size());
        if (finite_averaged || finite_averaged_friction) {
            struct ConstraintBuilder final {
                std::map<std::size_t, double> secondary;
                std::map<std::size_t, CartesianPoint3> coordinates;
                std::map<std::size_t, double> secondary_output;
                std::vector<AbaqusAveragedConstraint::FiniteSlidingSample> samples;
                CartesianPoint3 normal{};
                CartesianPoint3 tangent_first{};
                double area = 0.0;
            };

            std::vector<ConstraintBuilder> builders(secondary.boundary.nodes.size());
            const Matrix4& averaging = abaqus_quad4_averaging();
            for (std::size_t secondary_face_index = 0; secondary_face_index < secondary_faces.size();
                ++secondary_face_index) {
                const SecondaryContactFace& face = secondary_faces[secondary_face_index];
                const Quad4FaceElement& source_face = secondary.boundary.faces[secondary_face_index];
                for (std::size_t local_constraint = 0; local_constraint < 4; ++local_constraint) {
                    const auto found = std::find(secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(),
                        source_face.nodes[local_constraint]);
                    if (found == secondary.boundary.nodes.end())
                        throw std::logic_error("HEX8 finite-sliding averaged-contact node mapping failed");
                    ConstraintBuilder& builder =
                        builders[static_cast<std::size_t>(found - secondary.boundary.nodes.begin())];
                    const std::array<double, 2>& location = abaqus_quad4_constraint_locations()[local_constraint];
                    const Quad4FaceQuadraturePoint point =
                        make_quad4_face_quadrature_point(face.coordinates, location[0], location[1], 1.0);
                    const double local_area = point.weighted_measure;
                    builder.area += local_area;
                    for (std::size_t node = 0; node < 4; ++node) {
                        const double value = local_area * averaging[local_constraint * 4 + node];
                        builder.secondary[face.nodes[node]] += value;
                        const auto output_node = std::find(
                            secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(), source_face.nodes[node]);
                        if (output_node == secondary.boundary.nodes.end())
                            throw std::logic_error("HEX8 finite-sliding averaged-contact output mapping failed");
                        builder.secondary_output[static_cast<std::size_t>(
                            output_node - secondary.boundary.nodes.begin())] += value;
                        builder.coordinates[face.nodes[node]] = face.coordinates[node];
                    }
                    double minimum_distance = std::numeric_limits<double>::infinity();
                    std::size_t selected = std::numeric_limits<std::size_t>::max();
                    for (std::size_t primary_index = 0; primary_index < primary_faces.size(); ++primary_index) {
                        const PrimaryContactFace& primary_face = primary_faces[primary_index];
                        const Quad4ReferenceProjectionValue projection =
                            compute_quad4_reference_projection(face.coordinates, primary_face.coordinates, point.shape,
                                primary_material_orientation(primary_face.coordinates, primary_face.parent_centroid));
                        if (projection.projected && std::abs(projection.gap) < minimum_distance) {
                            minimum_distance = std::abs(projection.gap);
                            selected = primary_index;
                        }
                    }
                    if (selected == std::numeric_limits<std::size_t>::max())
                        throw std::invalid_argument("HEX8 finite-sliding surface contact '" + definition.name +
                                                    "' has an unprojected reference constraint center");
                    const double orientation = normal_orientation(
                        face.coordinates, face.parent_centroid, primary_faces[selected].parent_centroid);
                    const Quad4FaceQuadraturePoint normal_point = make_quad4_face_quadrature_point(
                        face.coordinates, (4.0 / 3.0) * location[0], (4.0 / 3.0) * location[1], 1.0);
                    CartesianPoint3 face_normal = cross(normal_point.tangent_xi, normal_point.tangent_eta);
                    const double normal_measure = std::sqrt(dot(face_normal, face_normal));
                    builder.normal.x += local_area * orientation * face_normal.x / normal_measure;
                    builder.normal.y += local_area * orientation * face_normal.y / normal_measure;
                    builder.normal.z += local_area * orientation * face_normal.z / normal_measure;
                    CartesianPoint3 tangent = point.tangent_xi;
                    const double tangent_measure = std::sqrt(dot(tangent, tangent));
                    if (!(tangent_measure > 0.0))
                        throw std::invalid_argument("HEX8 finite-sliding averaged contact has an undefined tangent");
                    const double tangent_orientation =
                        builder.area == local_area || dot(builder.tangent_first, tangent) >= 0.0 ? 1.0 : -1.0;
                    builder.tangent_first.x += local_area * tangent_orientation * tangent.x / tangent_measure;
                    builder.tangent_first.y += local_area * tangent_orientation * tangent.y / tangent_measure;
                    builder.tangent_first.z += local_area * tangent_orientation * tangent.z / tangent_measure;
                    builder.samples.push_back(
                        {secondary_face_index, local_constraint, orientation, tangent_orientation, {}, {}});
                }
            }
            std::vector<std::size_t> primary_nodes;
            std::map<std::size_t, CartesianPoint3> primary_coordinates;
            for (const PrimaryContactFace& face : primary_faces)
                for (std::size_t node = 0; node < 4; ++node) {
                    primary_nodes.push_back(face.nodes[node]);
                    primary_coordinates[face.nodes[node]] = face.coordinates[node];
                }
            std::sort(primary_nodes.begin(), primary_nodes.end());
            primary_nodes.erase(std::unique(primary_nodes.begin(), primary_nodes.end()), primary_nodes.end());
            for (std::size_t output = 0; output < builders.size(); ++output) {
                const ConstraintBuilder& builder = builders[output];
                if (!(builder.area > 0.0))
                    throw std::invalid_argument("HEX8 finite-sliding averaged contact has a nonpositive area");
                const double normal_measure = std::sqrt(dot(builder.normal, builder.normal));
                if (!(normal_measure > 0.0))
                    throw std::invalid_argument("HEX8 finite-sliding averaged contact has an undefined normal");
                AbaqusAveragedConstraint constraint{};
                constraint.contact = contact_value;
                constraint.secondary = output;
                constraint.history = output;
                constraint.finite_sliding = true;
                constraint.finite_region_normal = finite_region_normal;
                constraint.area = builder.area;
                constraint.normal = {builder.normal.x / normal_measure, builder.normal.y / normal_measure,
                    builder.normal.z / normal_measure};
                const double tangent_normal_component = dot(builder.tangent_first, constraint.normal);
                CartesianPoint3 tangent = {builder.tangent_first.x - tangent_normal_component * constraint.normal.x,
                    builder.tangent_first.y - tangent_normal_component * constraint.normal.y,
                    builder.tangent_first.z - tangent_normal_component * constraint.normal.z};
                const double tangent_measure = std::sqrt(dot(tangent, tangent));
                if (!(tangent_measure > 0.0))
                    throw std::invalid_argument(
                        "HEX8 finite-sliding averaged contact has an undefined averaged tangent");
                constraint.tangent_first = {
                    tangent.x / tangent_measure, tangent.y / tangent_measure, tangent.z / tangent_measure};
                constraint.reference_normal = constraint.normal;
                constraint.reference_tangent_first = constraint.tangent_first;
                constraint.finite_sliding_samples = builder.samples;
                for (const auto& entry : builder.secondary) {
                    constraint.nodes.push_back(entry.first);
                    constraint.gap_coefficients.push_back(-entry.second / builder.area);
                    constraint.reference_coordinates.push_back(builder.coordinates.at(entry.first));
                }
                for (const std::size_t node : primary_nodes) {
                    constraint.nodes.push_back(node);
                    constraint.gap_coefficients.push_back(0.0);
                    constraint.reference_coordinates.push_back(primary_coordinates.at(node));
                }
                for (const auto& entry : builder.secondary_output) {
                    constraint.secondary_output_nodes.push_back(entry.first);
                    constraint.secondary_coefficients.push_back(entry.second / builder.area);
                }
                constraint.tangent_first_coefficients.resize(constraint.nodes.size());
                constraint.tangent_second_coefficients.resize(constraint.nodes.size());
                constraint.traction_first_coefficients.resize(constraint.nodes.size());
                constraint.traction_second_coefficients.resize(constraint.nodes.size());
                constraint.secondary_tangent_first_coefficients.resize(constraint.secondary_output_nodes.size());
                constraint.secondary_tangent_second_coefficients.resize(constraint.secondary_output_nodes.size());
                CartesianPoint3 separation{};
                for (std::size_t node = 0; node < constraint.nodes.size(); ++node) {
                    const CartesianPoint3& point = constraint.reference_coordinates[node];
                    separation.x += constraint.gap_coefficients[node] * point.x;
                    separation.y += constraint.gap_coefficients[node] * point.y;
                    separation.z += constraint.gap_coefficients[node] * point.z;
                }
                constraint.reference_gap = dot(separation, constraint.normal);
                if (finite_averaged_friction) {
                    AbaqusAveragedConstraint friction = constraint;
                    friction.finite_region_normal = false;
                    friction.friction_only = true;
                    friction.history = _contact_histories[contact_value].size();
                    _contact_histories[contact_value].emplace_back();
                    _abaqus_averaged_constraints.push_back(std::move(constraint));
                    _abaqus_averaged_constraints.push_back(std::move(friction));
                } else {
                    _abaqus_averaged_constraints.push_back(std::move(constraint));
                }
            }
        }
        if (definition.mechanical && surface_to_surface && !finite_sliding) {
            struct ConstraintBuilder final {
                std::map<std::size_t, double> secondary, primary;
                std::map<std::size_t, CartesianPoint3> coordinates;
                std::map<std::size_t, double> secondary_output;
                CartesianPoint3 normal{};
                double area = 0.0;
            };

            struct Cell final {
                double xi_lower, xi_upper, eta_lower, eta_upper;
                std::size_t depth;
            };

            constexpr std::array<double, 3> points = {-0.7745966692414834, 0.0, 0.7745966692414834};
            constexpr std::array<double, 3> weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
            std::vector<ConstraintBuilder> builders(secondary.boundary.nodes.size());
            const Matrix4& averaging = abaqus_quad4_averaging();
            for (std::size_t secondary_face_index = 0; secondary_face_index < secondary_faces.size();
                ++secondary_face_index) {
                const SecondaryContactFace& face = secondary_faces[secondary_face_index];
                const Quad4FaceElement& source_face = secondary.boundary.faces[secondary_face_index];
                const auto projection = [&](double xi, double eta, std::size_t primary_index) {
                    const Quad4FaceQuadraturePoint point =
                        make_quad4_face_quadrature_point(face.coordinates, xi, eta, 1.0);
                    const PrimaryContactFace& primary_face = primary_faces.at(primary_index);
                    return compute_quad4_reference_projection(face.coordinates, primary_face.coordinates, point.shape,
                        normal_orientation(
                            primary_face.coordinates, face.parent_centroid, primary_face.parent_centroid));
                };
                const auto owner = [&](double xi, double eta) {
                    double minimum_distance = std::numeric_limits<double>::infinity();
                    std::size_t selected = std::numeric_limits<std::size_t>::max();
                    for (std::size_t primary_index = 0; primary_index < primary_faces.size(); ++primary_index) {
                        const Quad4ReferenceProjectionValue reference = projection(xi, eta, primary_index);
                        if (!reference.projected) continue;
                        const double distance = std::abs(reference.gap);
                        if (distance < minimum_distance || (distance == minimum_distance && primary_index < selected)) {
                            minimum_distance = distance;
                            selected = primary_index;
                        }
                    }
                    return selected;
                };

                Matrix4 gram{};
                double face_area = 0.0;
                for (const Quad4FaceQuadraturePoint& point : face.geometry.points) {
                    face_area += point.weighted_measure;
                    for (std::size_t row = 0; row < 4; ++row)
                        for (std::size_t column = 0; column < 4; ++column)
                            gram[row * 4 + column] += point.weighted_measure * point.shape[row] * point.shape[column];
                }
                if (!std::isfinite(face_area) || !(face_area > 0.0))
                    throw std::invalid_argument("HEX8 averaged contact has a nonpositive secondary face area");
                Matrix4 test_coefficients{};
                for (std::size_t local_constraint = 0; local_constraint < 4; ++local_constraint) {
                    Vector4 weighted_averaging{};
                    for (std::size_t node = 0; node < 4; ++node)
                        weighted_averaging[node] = 0.25 * face_area * averaging[local_constraint * 4 + node];
                    const Vector4 coefficients = solve4(gram, weighted_averaging);
                    for (std::size_t node = 0; node < 4; ++node)
                        test_coefficients[local_constraint * 4 + node] = coefficients[node];

                    const auto found = std::find(secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(),
                        source_face.nodes[local_constraint]);
                    if (found == secondary.boundary.nodes.end())
                        throw std::logic_error("HEX8 averaged contact secondary-node mapping failed");
                    ConstraintBuilder& builder =
                        builders[static_cast<std::size_t>(found - secondary.boundary.nodes.begin())];
                    const double local_area = 0.25 * face_area;
                    builder.area += local_area;
                    for (std::size_t node = 0; node < 4; ++node) {
                        const double value = weighted_averaging[node];
                        builder.secondary[face.nodes[node]] += value;
                        const auto output_node = std::find(
                            secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(), source_face.nodes[node]);
                        if (output_node == secondary.boundary.nodes.end())
                            throw std::logic_error("HEX8 averaged contact output-node mapping failed");
                        builder.secondary_output[static_cast<std::size_t>(
                            output_node - secondary.boundary.nodes.begin())] += value;
                        builder.coordinates[face.nodes[node]] = face.coordinates[node];
                    }
                    const std::array<double, 2>& location = abaqus_quad4_constraint_locations()[local_constraint];
                    const Quad4FaceQuadraturePoint constraint_point =
                        make_quad4_face_quadrature_point(face.coordinates, location[0], location[1], 1.0);
                    CartesianPoint3 face_normal = cross(constraint_point.tangent_xi, constraint_point.tangent_eta);
                    const std::size_t primary_index = owner(location[0], location[1]);
                    if (primary_index == std::numeric_limits<std::size_t>::max())
                        throw std::invalid_argument("HEX8 small-sliding surface contact '" + definition.name +
                                                    "' has an unprojected reference constraint center");
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
                        throw std::invalid_argument("HEX8 small-sliding surface contact '" + definition.name +
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
                            const std::size_t primary_index = owner(xi, eta);
                            if (primary_index == std::numeric_limits<std::size_t>::max())
                                throw std::invalid_argument("HEX8 small-sliding surface contact '" + definition.name +
                                                            "' has an unprojected reference integration point");
                            const Quad4FaceQuadraturePoint point =
                                make_quad4_face_quadrature_point(face.coordinates, xi, eta, weight);
                            const Quad4ReferenceProjectionValue reference = projection(xi, eta, primary_index);
                            if (!reference.projected)
                                throw std::logic_error(
                                    "HEX8 small-sliding contact lost its reference projection while building");
                            for (std::size_t local_constraint = 0; local_constraint < 4; ++local_constraint) {
                                double test = 0.0;
                                for (std::size_t node = 0; node < 4; ++node)
                                    test += test_coefficients[local_constraint * 4 + node] * point.shape[node];
                                const auto found = std::find(secondary.boundary.nodes.begin(),
                                    secondary.boundary.nodes.end(), source_face.nodes[local_constraint]);
                                ConstraintBuilder& builder =
                                    builders[static_cast<std::size_t>(found - secondary.boundary.nodes.begin())];
                                const PrimaryContactFace& primary_face = primary_faces[primary_index];
                                for (std::size_t node = 0; node < 4; ++node) {
                                    builder.primary[primary_face.nodes[node]] +=
                                        point.weighted_measure * test * reference.primary_shape[node];
                                    builder.coordinates[primary_face.nodes[node]] = primary_face.coordinates[node];
                                }
                            }
                        }
                }
            }
            for (std::size_t output = 0; output < builders.size(); ++output) {
                const ConstraintBuilder& builder = builders[output];
                if (!std::isfinite(builder.area) || !(builder.area > 0.0))
                    throw std::invalid_argument("HEX8 averaged contact has a nonpositive constraint area");
                const double normal_measure = std::sqrt(dot(builder.normal, builder.normal));
                if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
                    throw std::invalid_argument("HEX8 averaged contact has an undefined constraint normal");
                AbaqusAveragedConstraint constraint{};
                constraint.contact = contact_value;
                constraint.secondary = output;
                constraint.history = output;
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
                    throw std::invalid_argument("HEX8 averaged contact does not preserve rigid translation");
                CartesianPoint3 separation{};
                for (std::size_t node = 0; node < constraint.nodes.size(); ++node) {
                    const CartesianPoint3& coordinate = builder.coordinates.at(constraint.nodes[node]);
                    separation.x += constraint.gap_coefficients[node] * coordinate.x;
                    separation.y += constraint.gap_coefficients[node] * coordinate.y;
                    separation.z += constraint.gap_coefficients[node] * coordinate.z;
                }
                constraint.reference_gap = dot(separation, constraint.normal);
                _abaqus_averaged_constraints.push_back(std::move(constraint));
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
        if (definition.mechanical_sliding == MechanicalContactSliding::finite && !surface_to_surface)
            throw std::invalid_argument("finite sliding requires HEX20 surface_to_surface contact: " + definition.name);
        // Abaqus-style node-centered averaged constraints are used for
        // small-strain small sliding and, through a separate current-geometry
        // refresh path, finite sliding. Each curved constraint forms its current
        // normal and tangent from its own contributing face samples before the
        // node-centered friction update.
        const bool finite_sliding = definition.mechanical_sliding == MechanicalContactSliding::finite;
        const bool finite_averaged = surface_to_surface && finite_sliding,
                   abaqus_averaged =
                       finite_averaged ||
                       (surface_to_surface && !finite_sliding &&
                           _definition.regions[primary.region].strain_formulation == StrainFormulation::small &&
                           _definition.regions[secondary.region].strain_formulation == StrainFormulation::small);
        if (definition.friction_slip_tolerance > 0.0 && !surface_to_surface)
            throw std::invalid_argument("slip_tolerance requires HEX20 surface_to_surface contact: " + definition.name);
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
        _thermal_properties.push_back(
            {definition.thermal ? definition.gap_conductivity : 1.0, definition.thermal ? definition.minimum_gap : 1.0,
                definition.thermal ? definition.gap_heat_conductance_law : GapHeatConductanceLaw::gas_gap,
                definition.thermal ? definition.gap_conductance : 0.0,
                definition.thermal ? definition.gap_conductance_clearance_derivative : 0.0,
                definition.thermal ? definition.gap_conductance_pressure_derivative : 0.0,
                definition.thermal ? definition.gap_conductance_temperature_derivative : 0.0,
                definition.thermal ? definition.gap_conductance_reference_temperature : 0.0,
                definition.mechanical ? definition.penalty : 0.0});
        const double slip_tolerance =
            definition.mechanical && definition.friction_coefficient > 0.0 && surface_to_surface
                ? (definition.friction_slip_tolerance > 0.0 ? definition.friction_slip_tolerance : 5.0e-3)
                : 0.0;
        double reference_secondary_area = 0.0;
        if (slip_tolerance > 0.0)
            for (const Quad8FaceElement& face : secondary.boundary.faces)
                for (const Quad8FaceMechanicalQuadraturePoint& point :
                    make_quad8_face_geometry(face_coordinates(secondary_mesh, face)).mechanical_points)
                    reference_secondary_area += point.quadrature_weight * reference_measure(point);
        const double maximum_elastic_slip =
            slip_tolerance > 0.0 ? slip_tolerance * std::sqrt(reference_secondary_area /
                                                              static_cast<double>(secondary.boundary.faces.size()))
                                 : 0.0;
        if (!std::isfinite(maximum_elastic_slip) || maximum_elastic_slip < 0.0)
            throw std::invalid_argument(
                "HEX20 contact slip_tolerance gives an invalid elastic slip: " + definition.name);
        _mechanical_properties.push_back({definition.mechanical ? definition.penalty : 1.0,
            definition.mechanical ? definition.friction_coefficient : 0.0, false, maximum_elastic_slip});
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
        if (definition.mechanical && surface_to_surface && !finite_averaged) {
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
                std::vector<AbaqusAveragedConstraint::FiniteSlidingSample> samples;
                CartesianPoint3 normal{};
                CartesianPoint3 tangent_first{};
                double area = 0.0;
            };

            std::vector<ConstraintBuilder> builders(secondary.boundary.displacement_nodes.size());
            const Matrix8& averaging = abaqus_quad8_averaging();
            for (std::size_t secondary_face_index = 0; secondary_face_index < secondary_faces.size();
                ++secondary_face_index) {
                const Hex20SecondaryContactFace& face = secondary_faces[secondary_face_index];
                const Quad8FaceElement& source_face = secondary.boundary.faces[secondary_face_index];
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
                std::size_t representative_primary = std::numeric_limits<std::size_t>::max();
                if (finite_averaged) {
                    const Quad8FaceMechanicalQuadraturePoint face_center =
                        make_quad8_face_mechanical_point(face.coordinates, 0.0, 0.0, 1.0);
                    double representative_distance = std::numeric_limits<double>::infinity();
                    for (std::size_t primary_index = 0; primary_index < primary_faces.size(); ++primary_index) {
                        const Hex20PrimaryContactFace& primary_face = primary_faces[primary_index];
                        const Quad8ReferenceProjectionValue reference = compute_quad8_reference_projection(
                            face.coordinates, primary_face.coordinates, face_center.displacement_shape,
                            normal_orientation(
                                primary_face.coordinates, face.parent_centroid, primary_face.parent_centroid));
                        if (!reference.projected) continue;
                        const double distance = std::abs(reference.gap);
                        if (distance < representative_distance ||
                            (distance == representative_distance && primary_index < representative_primary)) {
                            representative_distance = distance;
                            representative_primary = primary_index;
                        }
                    }
                } else {
                    representative_primary = face.contact_primary_faces.at(face.contact_primary_faces.size() / 2);
                }
                if (representative_primary == std::numeric_limits<std::size_t>::max())
                    throw std::invalid_argument("HEX20 averaged surface contact '" + definition.name +
                                                "' has an unprojected reference face center");
                double face_area = 0.0;
                for (const Quad8FaceMechanicalQuadraturePoint& point : face.geometry.mechanical_points)
                    face_area += point.quadrature_weight * reference_measure(point);
                if (!std::isfinite(face_area) || !(face_area > 0.0))
                    throw std::invalid_argument("HEX20 averaged contact has a nonpositive secondary face area");
                for (std::size_t local_constraint = 0; local_constraint < 8; ++local_constraint) {
                    const double fraction = local_constraint < 4 ? 1.0 / 24.0 : 5.0 / 24.0;
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
                        const double value = local_area * averaging[local_constraint * 8 + node];
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
                    const std::size_t primary_index = representative_primary;
                    const CartesianPoint3 direction =
                        subtract(primary_faces[primary_index].parent_centroid, face.parent_centroid);
                    double orientation = 1.0;
                    if (dot(face_normal, direction) < 0.0) {
                        face_normal.x = -face_normal.x;
                        face_normal.y = -face_normal.y;
                        face_normal.z = -face_normal.z;
                        orientation = -1.0;
                    }
                    const double normal_measure = std::sqrt(dot(face_normal, face_normal));
                    if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
                        throw std::invalid_argument("HEX20 averaged contact has an undefined local normal");
                    builder.normal.x += local_area * face_normal.x / normal_measure;
                    builder.normal.y += local_area * face_normal.y / normal_measure;
                    builder.normal.z += local_area * face_normal.z / normal_measure;
                    CartesianPoint3 tangent = constraint_point.tangent_xi;
                    const double tangent_measure = std::sqrt(dot(tangent, tangent));
                    if (!std::isfinite(tangent_measure) || !(tangent_measure > 0.0))
                        throw std::invalid_argument("HEX20 averaged contact has an undefined local tangent");
                    const double tangent_orientation =
                        builder.samples.empty() || dot(builder.tangent_first, tangent) >= 0.0 ? 1.0 : -1.0;
                    builder.tangent_first.x += local_area * tangent_orientation * tangent.x / tangent_measure;
                    builder.tangent_first.y += local_area * tangent_orientation * tangent.y / tangent_measure;
                    builder.tangent_first.z += local_area * tangent_orientation * tangent.z / tangent_measure;
                    if (finite_averaged) {
                        builder.samples.push_back(
                            {secondary_face_index, local_constraint, orientation, tangent_orientation, {}, {}});
                    } else {
                        for (const AbaqusQuad8TransferSample& sample :
                            abaqus_quad8_primary_transfer_rule(local_constraint)) {
                            const double xi = 2.0 * sample.first - 1.0, eta = 2.0 * sample.second - 1.0;
                            const std::size_t primary_index = owner(xi, eta);
                            if (primary_index == std::numeric_limits<std::size_t>::max())
                                throw std::invalid_argument("HEX20 small-sliding surface contact '" + definition.name +
                                                            "' has an unprojected Abaqus primary-transfer sample");
                            const Quad8ReferenceProjectionValue reference = projection(xi, eta, primary_index);
                            if (!reference.projected)
                                throw std::logic_error(
                                    "HEX20 small-sliding contact lost an Abaqus primary-transfer projection");
                            const Hex20PrimaryContactFace& primary_face = primary_faces[primary_index];
                            for (std::size_t node = 0; node < 8; ++node) {
                                builder.primary[primary_face.displacement_nodes[node]] +=
                                    local_area * sample.weight * reference.primary_shape[node];
                                builder.coordinates[primary_face.displacement_nodes[node]] =
                                    primary_face.coordinates[node];
                            }
                        }
                    }
                }
            }
            std::vector<std::size_t> finite_primary_nodes;
            std::map<std::size_t, CartesianPoint3> finite_primary_coordinates;
            if (finite_averaged) {
                for (const Hex20PrimaryContactFace& face : primary_faces)
                    for (std::size_t node = 0; node < 8; ++node) {
                        finite_primary_nodes.push_back(face.displacement_nodes[node]);
                        finite_primary_coordinates[face.displacement_nodes[node]] = face.coordinates[node];
                    }
                std::sort(finite_primary_nodes.begin(), finite_primary_nodes.end());
                finite_primary_nodes.erase(
                    std::unique(finite_primary_nodes.begin(), finite_primary_nodes.end()), finite_primary_nodes.end());
            }
            for (std::size_t output = 0; output < builders.size(); ++output) {
                const ConstraintBuilder& builder = builders[output];
                if (!std::isfinite(builder.area) || !(builder.area > 0.0))
                    throw std::invalid_argument("HEX20 averaged contact has a nonpositive constraint area");
                const double normal_measure = std::sqrt(dot(builder.normal, builder.normal));
                if (!std::isfinite(normal_measure) || !(normal_measure > 0.0))
                    throw std::invalid_argument("HEX20 averaged contact has an undefined constraint normal");
                AbaqusAveragedConstraint constraint{};
                constraint.contact = contact_value;
                constraint.secondary = output;
                constraint.history = output;
                constraint.normal = {builder.normal.x / normal_measure, builder.normal.y / normal_measure,
                    builder.normal.z / normal_measure};
                constraint.area = builder.area;
                if (finite_averaged) {
                    const double tangent_normal_component = dot(builder.tangent_first, constraint.normal);
                    CartesianPoint3 tangent = {builder.tangent_first.x - tangent_normal_component * constraint.normal.x,
                        builder.tangent_first.y - tangent_normal_component * constraint.normal.y,
                        builder.tangent_first.z - tangent_normal_component * constraint.normal.z};
                    const double tangent_measure = std::sqrt(dot(tangent, tangent));
                    if (!std::isfinite(tangent_measure) || !(tangent_measure > 0.0))
                        throw std::invalid_argument(
                            "HEX20 finite-sliding averaged contact has an undefined averaged tangent");
                    constraint.tangent_first = {
                        tangent.x / tangent_measure, tangent.y / tangent_measure, tangent.z / tangent_measure};
                    constraint.reference_normal = constraint.normal;
                    constraint.reference_tangent_first = constraint.tangent_first;
                    constraint.finite_sliding = true;
                    constraint.finite_sliding_samples = builder.samples;
                }
                double secondary_sum = 0.0, primary_sum = 0.0;
                for (const auto& entry : builder.secondary) {
                    const double coefficient = -entry.second / builder.area;
                    if (!finite_averaged && std::abs(coefficient) <= 1.0e-14) continue;
                    constraint.nodes.push_back(entry.first);
                    constraint.gap_coefficients.push_back(coefficient);
                    if (finite_averaged)
                        constraint.reference_coordinates.push_back(builder.coordinates.at(entry.first));
                    secondary_sum -= coefficient;
                }
                if (finite_averaged) {
                    for (const std::size_t node : finite_primary_nodes) {
                        constraint.nodes.push_back(node);
                        constraint.gap_coefficients.push_back(0.0);
                        constraint.reference_coordinates.push_back(finite_primary_coordinates.at(node));
                    }
                } else {
                    for (const auto& entry : builder.primary) {
                        const double coefficient = entry.second / builder.area;
                        if (std::abs(coefficient) <= 1.0e-14) continue;
                        constraint.nodes.push_back(entry.first);
                        constraint.gap_coefficients.push_back(coefficient);
                        primary_sum += coefficient;
                    }
                }
                for (const auto& entry : builder.secondary_output) {
                    constraint.secondary_output_nodes.push_back(entry.first);
                    constraint.secondary_coefficients.push_back(entry.second / builder.area);
                }
                if (finite_averaged) {
                    constraint.normal_gap_coefficients.resize(constraint.nodes.size());
                    constraint.secondary_normal_coefficients.resize(constraint.secondary_output_nodes.size());
                }
                if (std::abs(secondary_sum - 1.0) > 1.0e-10 ||
                    (!finite_averaged && std::abs(primary_sum - 1.0) > 1.0e-10))
                    throw std::invalid_argument("HEX20 averaged contact does not preserve rigid translation");
                CartesianPoint3 separation{};
                for (std::size_t node = 0; node < constraint.nodes.size(); ++node) {
                    const CartesianPoint3& coordinate = finite_averaged
                                                            ? constraint.reference_coordinates[node]
                                                            : builder.coordinates.at(constraint.nodes[node]);
                    separation.x += constraint.gap_coefficients[node] * coordinate.x;
                    separation.y += constraint.gap_coefficients[node] * coordinate.y;
                    separation.z += constraint.gap_coefficients[node] * coordinate.z;
                }
                constraint.reference_gap = dot(separation, constraint.normal);
                _abaqus_averaged_constraints.push_back(std::move(constraint));
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
    const bool surface_to_surface = _definition.contacts[metadata.contact].mechanical_discretization ==
                                    MechanicalContactDiscretization::surface_to_surface;
    MechanicalCandidate result{};
    result.contact = metadata.contact;
    result.nodes = nodes;
    result.surface_to_surface = surface_to_surface;
    result.secondary = metadata.secondary;
    result.primary = primary;
    result.secondary_face = metadata.secondary_face;
    result.secondary_local_point = metadata.secondary_local_point;
    if (surface_to_surface) {
        const Quad4FaceQuadraturePoint& contact_point = secondary.contact_points.at(metadata.secondary_local_point);
        const std::array<double, 2>& constraint_location =
            abaqus_quad4_constraint_locations().at(metadata.secondary_local_point);
        const Quad4FaceQuadraturePoint normal_point = make_quad4_face_quadrature_point(
            secondary.coordinates, (4.0 / 3.0) * constraint_location[0], (4.0 / 3.0) * constraint_location[1], 1.0);
        result.surface_geometry = {secondary.coordinates, primary_face.coordinates, contact_point.shape,
            contact_point.derivative_xi, contact_point.derivative_eta, normal_point.derivative_xi,
            normal_point.derivative_eta, 1.0,
            primary_material_orientation(primary_face.coordinates, primary_face.parent_centroid),
            normal_orientation(secondary.coordinates, secondary.parent_centroid, primary_face.parent_centroid)};
        return result;
    }
    std::array<std::array<double, 4>, 4> shapes{}, derivatives_xi{}, derivatives_eta{};
    for (std::size_t q = 0; q < secondary.geometry.points.size(); ++q) {
        shapes[q] = secondary.geometry.points[q].shape;
        derivatives_xi[q] = secondary.geometry.points[q].derivative_xi;
        derivatives_eta[q] = secondary.geometry.points[q].derivative_eta;
    }
    result.node_geometry = {secondary.coordinates, primary_face.coordinates, shapes, derivatives_xi, derivatives_eta,
        metadata.secondary_local_point,
        normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid)};
    return result;
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
    const bool finite_sliding =
        _definition.contacts[metadata.contact].mechanical_sliding == MechanicalContactSliding::finite;
    const double orientation =
        finite_sliding
            ? primary_material_orientation(primary_face.coordinates, primary_face.parent_centroid)
            : normal_orientation(primary_face.coordinates, secondary.parent_centroid, primary_face.parent_centroid);
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
            finite_sliding ? orientation : secondary.contact_normal_orientations.at(metadata.secondary_local_point),
            finite_sliding};
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
        bool thermal_touched = false, mechanical_touched = false;
        for (std::size_t point = _thermal_contact_offsets[contact]; point < _thermal_contact_offsets[contact + 1];
            ++point)
            thermal_touched = thermal_touched || _touched_thermal_points[point] != 0U;
        for (std::size_t node = _mechanical_contact_offsets[contact]; node < _mechanical_contact_offsets[contact + 1];
            ++node)
            mechanical_touched = mechanical_touched || _touched_mechanical_nodes[node] != 0U;
        if (!thermal_touched && !mechanical_touched) continue;
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
    const std::size_t point_count = _uses_hex20 ? _hex20_mechanical_points.size() : _mechanical_points.size();
    for (std::size_t entry = begin; entry < end; ++entry) {
        const std::size_t local = entry - ranges.mechanical_begin;
        if (local >= point_count) {
            const AbaqusAveragedConstraint& constraint = _abaqus_averaged_constraints.at(local - point_count);
            _touched_mechanical_nodes[mechanical_node_index(constraint.contact, constraint.history)] = 1U;
        } else if (_uses_hex20) {
            const Hex20MechanicalPoint& point = _hex20_mechanical_points[local];
            _touched_mechanical_nodes[mechanical_node_index(point.contact, point.secondary)] = 1U;
        } else {
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
            const bool finite_sliding =
                surface_to_surface &&
                _definition.contacts[metadata.contact].mechanical_sliding == MechanicalContactSliding::finite;
            if (surface_to_surface && !finite_sliding) {
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
        for (const AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
            const std::size_t node = mechanical_node_index(constraint.contact, constraint.history);
            if (_touched_mechanical_nodes[node] != 0U && constraint.projected)
                _mechanical_selected_primary[node] = constraint.primary_face;
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
            const Quad4SurfaceContactLocalValues candidate_state = contact_state(candidate.nodes, state);
            const ContactProjectionValue value =
                candidate.surface_to_surface
                    ? compute_quad4_to_quad4_contact_projection(candidate.surface_geometry, candidate_state)
                    : compute_node_to_quad4_contact_projection(candidate.node_geometry, candidate_state);
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
                representative.surface_to_surface
                    ? mechanical_search_point(representative.surface_geometry, representative_state)
                    : mechanical_search_point(representative.node_geometry, representative_state),
                _contact_search_query);
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
        const Quad4SurfaceContactLocalValues candidate_state = contact_state(candidate.nodes, state);
        const ContactProjectionValue value =
            candidate.surface_to_surface
                ? compute_quad4_to_quad4_contact_projection(candidate.surface_geometry, candidate_state)
                : compute_node_to_quad4_contact_projection(candidate.node_geometry, candidate_state);
        if (value.projected) _mechanical_active_primary[point] = primary;
    }
    for (const AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        const std::size_t node = mechanical_node_index(constraint.contact, constraint.history);
        if (_touched_mechanical_nodes[node] != 0U && constraint.projected)
            _mechanical_selected_primary[node] = constraint.primary_face;
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
            const ContributionRanges ranges = contribution_ranges();
            const std::size_t mechanical_begin = std::max(first, ranges.mechanical_begin),
                              mechanical_end = std::min(last, ranges.boundary_begin),
                              point_count = _hex20_mechanical_points.size();
            for (std::size_t entry = mechanical_begin; entry < mechanical_end; ++entry) {
                const std::size_t local = entry - ranges.mechanical_begin;
                if (local < point_count) continue;
                const AbaqusAveragedConstraint& constraint = _abaqus_averaged_constraints.at(local - point_count);
                if (!constraint.finite_sliding) continue;
                touched_contacts[constraint.contact] = 1U;
                for (const AbaqusAveragedConstraint::FiniteSlidingSample& sample : constraint.finite_sliding_samples) {
                    const Hex20SecondaryContactFace& face =
                        _hex20_secondary_contact_faces[constraint.contact][sample.secondary_face];
                    append_hex20_face(face.temperature_nodes, face.displacement_nodes);
                }
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
    _fully_validated_contact_state_current = false;
    mark_touched_thermal_points(first, last);
    mark_touched_mechanical_nodes(first, last);
    update_contact_search_trees(state);
    update_thermal_candidates(first, last, state);
    refresh_finite_averaged_constraints(state);
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
            const bool natural_finite_surface_release =
                !_uses_hex20 &&
                _definition.contacts[contact].mechanical_discretization ==
                    MechanicalContactDiscretization::surface_to_surface &&
                _definition.contacts[contact].mechanical_sliding == MechanicalContactSliding::finite &&
                _definition.contacts[contact].friction_coefficient == 0.0 &&
                std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
                    [contact](const AbaqusAveragedConstraint& constraint) {
                        return constraint.contact == contact && constraint.finite_sliding && !constraint.friction_only;
                    });
            if (unprojected != 0 && !natural_finite_surface_release)
                throw std::domain_error(
                    "Three-dimensional mechanical contact '" + _definition.contacts[contact].name +
                    "' lost projection for " + std::to_string(unprojected) +
                    (_definition.contacts[contact].mechanical_discretization ==
                                MechanicalContactDiscretization::surface_to_surface
                            ? (_definition.contacts[contact].mechanical_sliding == MechanicalContactSliding::finite
                                      ? " finite-sliding integration points after searching every primary face"
                                      : (_uses_hex20 ? " fixed small-sliding integration points"
                                                     : " fixed small-sliding averaged constraints"))
                            : " secondary nodes after searching every primary face"));
        }
    }
    if (first == 0 && last == contribution_count()) {
        _fully_validated_contact_state = state;
        _fully_validated_contact_state_current = true;
    }
}

double SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    if (state.size() != dof_count())
        throw std::invalid_argument("Three-dimensional committed contact state size mismatch");
    const bool validated_state_matches = _fully_validated_contact_state_current &&
                                         _fully_validated_contact_state.size() == state.size() &&
                                         std::equal(state.begin(), state.end(), _fully_validated_contact_state.begin());
    if (!validated_state_matches) {
        update_contact_search_trees(state);
        refresh_finite_averaged_constraints(state);
        update_mechanical_candidates(0, contribution_count(), state);
        _fully_validated_contact_state = state;
        _fully_validated_contact_state_current = true;
    }
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    double friction_dissipation = 0.0;
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
            trial.cartesian_total_tangential_slip = value.tangential_slip;
            trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
            if (candidate.surface_to_surface && _mechanical_properties[candidate.contact].friction_coefficient > 0.0) {
                trial.cartesian_tangent_basis_initialized = true;
                trial.cartesian_contact_normal = value.normal;
                trial.cartesian_contact_tangent_first = value.tangent_first;
            }
            if (updated[candidate.contact][candidate.secondary]) {
                const ContactPointHistory& prior = staged[candidate.contact][candidate.secondary];
                bool equal = prior.sliding == trial.sliding;
                for (std::size_t component = 0; component < 3; ++component) {
                    const double scale = std::max({1.0, std::abs(prior.cartesian_elastic_tangential_slip[component]),
                        std::abs(trial.cartesian_elastic_tangential_slip[component])});
                    equal = equal && std::abs(prior.cartesian_elastic_tangential_slip[component] -
                                              trial.cartesian_elastic_tangential_slip[component]) <=
                                         64.0 * std::numeric_limits<double>::epsilon() * scale;
                    const double total_scale =
                        std::max({1.0, std::abs(prior.cartesian_total_tangential_slip[component]),
                            std::abs(trial.cartesian_total_tangential_slip[component])});
                    equal = equal && std::abs(prior.cartesian_total_tangential_slip[component] -
                                              trial.cartesian_total_tangential_slip[component]) <=
                                         64.0 * std::numeric_limits<double>::epsilon() * total_scale;
                    equal = equal &&
                            prior.cartesian_contact_normal[component] == trial.cartesian_contact_normal[component] &&
                            prior.cartesian_contact_tangent_first[component] ==
                                trial.cartesian_contact_tangent_first[component];
                }
                equal = equal && prior.cartesian_tangent_basis_initialized == trial.cartesian_tangent_basis_initialized;
                if (!equal) throw std::logic_error("HEX20 secondary faces disagree on contact-node friction history");
                continue;
            }
            staged[candidate.contact][candidate.secondary] = trial;
            updated[candidate.contact][candidate.secondary] = true;
            friction_dissipation += value.friction_dissipation;
        }
        for (const AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
            std::vector<std::size_t> dofs;
            averaged_constraint_dofs(constraint, dofs);
            std::vector<double> current(dofs.size()), committed(dofs.size());
            for (std::size_t local = 0; local < dofs.size(); ++local) {
                current[local] = state.at(dofs[local]);
                committed[local] = _committed_contact_solution.at(dofs[local]);
            }
            const AbaqusAveragedConstraintValue value = averaged_constraint_value(
                constraint, current, committed, _contact_histories[constraint.contact][constraint.history]);
            ContactPointHistory trial = _contact_histories[constraint.contact][constraint.history];
            trial.sliding = value.sliding;
            trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
            trial.cartesian_total_tangential_slip = value.tangential_slip;
            if (constraint.finite_sliding && _mechanical_properties[constraint.contact].friction_coefficient > 0.0) {
                trial.cartesian_tangent_basis_initialized = true;
                trial.cartesian_contact_normal =
                    std::array<double, 3>{constraint.normal.x, constraint.normal.y, constraint.normal.z};
                trial.cartesian_contact_tangent_first = value.tangent_first;
            }
            if (updated[constraint.contact][constraint.history])
                throw std::logic_error("Abaqus-style averaged constraints share one friction-history slot");
            staged[constraint.contact][constraint.history] = trial;
            updated[constraint.contact][constraint.history] = true;
            friction_dissipation += value.friction_dissipation;
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
        return friction_dissipation;
    }
    for (std::size_t point = 0; point < _mechanical_active_primary.size(); ++point) {
        const std::size_t primary = _mechanical_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        const Quad4SurfaceContactLocalValues current = contact_state(candidate.nodes, state),
                                             committed = contact_state(candidate.nodes, _committed_contact_solution);
        NormalContactProperties point_properties = _mechanical_properties[candidate.contact];
        if (candidate.surface_to_surface &&
            std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
                [&candidate](const AbaqusAveragedConstraint& constraint) {
                    return constraint.contact == candidate.contact && constraint.friction_only;
                })) {
            point_properties.friction_coefficient = 0.0;
            point_properties.maximum_elastic_slip = 0.0;
        }
        const CartesianContactPointValue value =
            candidate.surface_to_surface
                ? compute_quad4_to_quad4_contact_value(point_properties, candidate.surface_geometry, current, committed,
                      _contact_histories[candidate.contact][candidate.secondary])
                : compute_node_to_quad4_contact_value(_mechanical_properties[candidate.contact],
                      candidate.node_geometry, current, committed,
                      _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        ContactPointHistory trial = _contact_histories[candidate.contact][candidate.secondary];
        trial.sliding = value.sliding;
        trial.cartesian_total_tangential_slip = value.tangential_slip;
        trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
        if (candidate.surface_to_surface && point_properties.friction_coefficient > 0.0) {
            trial.cartesian_tangent_basis_initialized = true;
            trial.cartesian_contact_normal = value.normal;
            trial.cartesian_contact_tangent_first = value.tangent_first;
        }
        if (updated[candidate.contact][candidate.secondary]) {
            const ContactPointHistory& prior = staged[candidate.contact][candidate.secondary];
            bool equal = prior.sliding == trial.sliding;
            for (std::size_t component = 0; component < 3; ++component) {
                const double scale = std::max({1.0, std::abs(prior.cartesian_elastic_tangential_slip[component]),
                    std::abs(trial.cartesian_elastic_tangential_slip[component])});
                equal = equal && std::abs(prior.cartesian_elastic_tangential_slip[component] -
                                          trial.cartesian_elastic_tangential_slip[component]) <=
                                     64.0 * std::numeric_limits<double>::epsilon() * scale;
                const double total_scale = std::max({1.0, std::abs(prior.cartesian_total_tangential_slip[component]),
                    std::abs(trial.cartesian_total_tangential_slip[component])});
                equal = equal && std::abs(prior.cartesian_total_tangential_slip[component] -
                                          trial.cartesian_total_tangential_slip[component]) <=
                                     64.0 * std::numeric_limits<double>::epsilon() * total_scale;
                equal = equal &&
                        prior.cartesian_contact_normal[component] == trial.cartesian_contact_normal[component] &&
                        prior.cartesian_contact_tangent_first[component] ==
                            trial.cartesian_contact_tangent_first[component];
            }
            equal = equal && prior.cartesian_tangent_basis_initialized == trial.cartesian_tangent_basis_initialized;
            if (!equal)
                throw std::logic_error(
                    "Three-dimensional secondary half-faces disagree on contact-node friction history");
            continue;
        }
        staged[candidate.contact][candidate.secondary] = trial;
        updated[candidate.contact][candidate.secondary] = true;
        friction_dissipation += value.friction_dissipation;
    }
    for (const AbaqusAveragedConstraint& constraint : _abaqus_averaged_constraints) {
        std::vector<std::size_t> dofs;
        averaged_constraint_dofs(constraint, dofs);
        std::vector<double> current(dofs.size()), committed(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local) {
            current[local] = state.at(dofs[local]);
            committed[local] = _committed_contact_solution.at(dofs[local]);
        }
        const AbaqusAveragedConstraintValue value = averaged_constraint_value(
            constraint, current, committed, _contact_histories[constraint.contact][constraint.history]);
        ContactPointHistory trial = _contact_histories[constraint.contact][constraint.history];
        trial.sliding = value.sliding;
        trial.cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
        trial.cartesian_total_tangential_slip = value.tangential_slip;
        if (constraint.finite_sliding && _mechanical_properties[constraint.contact].friction_coefficient > 0.0) {
            trial.cartesian_tangent_basis_initialized = true;
            trial.cartesian_contact_normal =
                std::array<double, 3>{constraint.normal.x, constraint.normal.y, constraint.normal.z};
            trial.cartesian_contact_tangent_first = value.tangent_first;
        }
        if (updated[constraint.contact][constraint.history])
            throw std::logic_error("Abaqus-style averaged constraints share one friction-history slot");
        staged[constraint.contact][constraint.history] = trial;
        updated[constraint.contact][constraint.history] = true;
        friction_dissipation += value.friction_dissipation;
    }
    for (std::size_t contact = 0; contact < _definition.contacts.size(); ++contact) {
        if (!_definition.contacts[contact].mechanical) continue;
        if (std::find(updated[contact].begin(), updated[contact].end(), false) != updated[contact].end())
            throw std::domain_error(
                _definition.contacts[contact].mechanical_discretization ==
                        MechanicalContactDiscretization::surface_to_surface
                    ? "Cannot commit HEX8 friction history for an unprojected finite-sliding contact point"
                    : "Cannot commit three-dimensional friction history for an unprojected node");
    }
    _contact_histories.swap(staged);
    _committed_contact_solution = state;
    return friction_dissipation;
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
            for (double component : history.cartesian_total_tangential_slip)
                if (!std::isfinite(component))
                    throw std::invalid_argument("Three-dimensional restored total-slip history is invalid");
            for (double component : history.cartesian_contact_normal)
                if (!std::isfinite(component))
                    throw std::invalid_argument("Three-dimensional restored contact normal is invalid");
            for (double component : history.cartesian_contact_tangent_first)
                if (!std::isfinite(component))
                    throw std::invalid_argument("Three-dimensional restored contact tangent is invalid");
            if (history.cartesian_tangent_basis_initialized) {
                double normal_norm = 0.0, tangent_norm = 0.0, orthogonality = 0.0;
                for (std::size_t component = 0; component < 3; ++component) {
                    normal_norm +=
                        history.cartesian_contact_normal[component] * history.cartesian_contact_normal[component];
                    tangent_norm += history.cartesian_contact_tangent_first[component] *
                                    history.cartesian_contact_tangent_first[component];
                    orthogonality += history.cartesian_contact_normal[component] *
                                     history.cartesian_contact_tangent_first[component];
                }
                constexpr double basis_tolerance = 1.0e-8;
                if (std::abs(normal_norm - 1.0) > basis_tolerance || std::abs(tangent_norm - 1.0) > basis_tolerance ||
                    std::abs(orthogonality) > basis_tolerance)
                    throw std::invalid_argument("Three-dimensional restored contact basis is not orthonormal");
            }
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
                std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, {}, {}, {}, false});
        }
        if (summarize_averaged_contact(contact_value, state, result)) return result;
        if (_definition.contacts[contact_value].mechanical_discretization ==
            MechanicalContactDiscretization::surface_to_surface) {
            std::vector<std::array<double, 3>> weighted_total_slip(result.size()), weighted_elastic_slip(result.size());
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
                    for (std::size_t component = 0; component < 3; ++component) {
                        weighted_total_slip[output_node][component] += nodal_area * value.tangential_slip[component];
                        weighted_elastic_slip[output_node][component] +=
                            nodal_area * value.elastic_tangential_slip[component];
                    }
                    summary.sliding = summary.sliding || value.sliding;
                }
            }
            NormalContactProperties recovery_properties = _mechanical_properties[contact_value];
            recovery_properties.friction_coefficient = 0.0;
            const bool finite_sliding =
                _definition.contacts[contact_value].mechanical_sliding == MechanicalContactSliding::finite;
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
                        if (finite_sliding) {
                            Hex20MechanicalCandidate trial{};
                            trial.contact = contact_value;
                            trial.secondary_temperature_nodes = secondary_geometry.temperature_nodes;
                            trial.primary_temperature_nodes = primary_geometry.temperature_nodes;
                            trial.secondary_displacement_nodes = secondary_geometry.displacement_nodes;
                            trial.primary_displacement_nodes = primary_geometry.displacement_nodes;
                            trial.surface_to_surface = true;
                            trial.surface_geometry = {secondary_geometry.coordinates, primary_geometry.coordinates,
                                quadrature.displacement_shape, quadrature.derivative_xi, quadrature.derivative_eta, {},
                                {}, {}, quadrature.quadrature_weight,
                                primary_material_orientation(
                                    primary_geometry.coordinates, primary_geometry.parent_centroid),
                                true};
                            const ContactProjectionValue projection = compute_quad8_to_quad8_contact_projection(
                                trial.surface_geometry, hex20_contact_state(trial, state));
                            if (!projection.projected || std::abs(projection.gap) >= minimum_distance) continue;
                            minimum_distance = std::abs(projection.gap);
                            selected_primary = primary_face;
                            continue;
                        }
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
                        finite_sliding ? primary_material_orientation(
                                             primary_geometry.coordinates, primary_geometry.parent_centroid)
                                       : normal_orientation(primary_geometry.coordinates,
                                             secondary_geometry.parent_centroid, primary_geometry.parent_centroid),
                        finite_sliding};
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
                summary.constraint_pressure = summary.pressure;
                summary.tangential_traction = summary.tangential_force / summary.tributary_area;
                for (std::size_t component = 0; component < 3; ++component) {
                    summary.tangential_slip[component] = weighted_total_slip[node][component] / summary.tributary_area;
                    summary.elastic_tangential_slip[component] =
                        weighted_elastic_slip[node][component] / summary.tributary_area;
                }
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
            summary.tangential_slip = value.tangential_slip;
            summary.sliding = value.sliding;
        }
        for (CartesianContactNodeSummary& summary : result)
            if (summary.tributary_area != 0.0) {
                summary.pressure = summary.contact_force / summary.tributary_area;
                summary.constraint_pressure = summary.pressure;
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
            std::numeric_limits<double>::infinity(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, {}, {}, {}, {}, false});
    }
    if (summarize_averaged_contact(contact_value, state, result)) return result;
    const bool separately_averaged_friction = std::any_of(_abaqus_averaged_constraints.begin(),
        _abaqus_averaged_constraints.end(), [contact_value](const AbaqusAveragedConstraint& constraint) {
            return constraint.contact == contact_value && constraint.friction_only;
        });
    std::vector<std::array<double, 3>> recovered_total_slip(result.size()), recovered_elastic_slip(result.size());
    std::vector<double> recovered_slip_weight(result.size());
    std::vector<std::array<std::array<double, 3>, 4>> contact_point_total_slip(secondary.boundary.faces.size()),
        contact_point_elastic_slip(secondary.boundary.faces.size());
    std::vector<std::array<unsigned char, 4>> contact_point_slip_set(secondary.boundary.faces.size());
    std::vector<double> contact_face_area(secondary.boundary.faces.size());
    for (std::size_t point = 0; point < _mechanical_active_primary.size(); ++point) {
        const std::size_t primary = _mechanical_active_primary[point];
        if (primary == std::numeric_limits<std::size_t>::max()) continue;
        const MechanicalCandidate candidate = mechanical_candidate(point, primary);
        if (candidate.contact != contact_value) continue;
        const std::size_t node = mechanical_node_index(candidate.contact, candidate.secondary);
        if (candidate.primary != _mechanical_selected_primary[node]) continue;
        const Quad4SurfaceContactLocalValues current = contact_state(candidate.nodes, state),
                                             committed = contact_state(candidate.nodes, _committed_contact_solution);
        NormalContactProperties point_properties = _mechanical_properties[candidate.contact];
        if (candidate.surface_to_surface &&
            std::any_of(_abaqus_averaged_constraints.begin(), _abaqus_averaged_constraints.end(),
                [&candidate](const AbaqusAveragedConstraint& constraint) {
                    return constraint.contact == candidate.contact && constraint.friction_only;
                })) {
            point_properties.friction_coefficient = 0.0;
            point_properties.maximum_elastic_slip = 0.0;
        }
        const CartesianContactPointValue value =
            candidate.surface_to_surface
                ? compute_quad4_to_quad4_contact_value(point_properties, candidate.surface_geometry, current, committed,
                      _contact_histories[candidate.contact][candidate.secondary])
                : compute_node_to_quad4_contact_value(_mechanical_properties[candidate.contact],
                      candidate.node_geometry, current, committed,
                      _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected) continue;
        if (candidate.surface_to_surface) {
            if (!separately_averaged_friction) {
                contact_point_total_slip[candidate.secondary_face][candidate.secondary_local_point] =
                    value.tangential_slip;
                contact_point_elastic_slip[candidate.secondary_face][candidate.secondary_local_point] =
                    value.elastic_tangential_slip;
                contact_point_slip_set[candidate.secondary_face][candidate.secondary_local_point] = 1U;
                contact_face_area[candidate.secondary_face] += value.tributary_area;
            }
            const Quad4FaceElement& face = secondary.boundary.faces.at(candidate.secondary_face);
            for (std::size_t local_node = 0; local_node < face.nodes.size(); ++local_node) {
                const auto found =
                    std::find(secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(), face.nodes[local_node]);
                if (found == secondary.boundary.nodes.end())
                    throw std::logic_error("HEX8 finite-sliding surface-contact output node mapping failed");
                const std::size_t output_node = static_cast<std::size_t>(found - secondary.boundary.nodes.begin());
                CartesianContactNodeSummary& summary = result[output_node];
                const double shape = candidate.surface_geometry.secondary_shape[local_node],
                             nodal_area = shape * value.tributary_area;
                summary.projected = true;
                summary.primary_face = std::min(summary.primary_face, candidate.primary);
                summary.gap = std::min(summary.gap, value.gap);
                summary.tributary_area += nodal_area;
                summary.contact_force += shape * value.contact_force;
                for (std::size_t component = 0; component < 3; ++component) {
                    summary.normal_contact_force[component] += shape * value.contact_force * value.normal[component];
                    summary.tangential_contact_force[component] +=
                        shape * value.tributary_area * value.tangential_traction_vector[component];
                }
                summary.tangential_force += shape * value.tangential_force;
                summary.sliding = summary.sliding || value.sliding;
            }
            continue;
        }
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
        summary.tangential_slip = value.tangential_slip;
        summary.sliding = value.sliding;
    }
    if (!separately_averaged_friction) {
        for (std::size_t face_index = 0; face_index < secondary.boundary.faces.size(); ++face_index) {
            if (std::find(contact_point_slip_set[face_index].begin(), contact_point_slip_set[face_index].end(), 0U) !=
                contact_point_slip_set[face_index].end())
                continue;
            const Quad4FaceElement& face = secondary.boundary.faces[face_index];
            for (std::size_t local_node = 0; local_node < 4; ++local_node) {
                const auto found =
                    std::find(secondary.boundary.nodes.begin(), secondary.boundary.nodes.end(), face.nodes[local_node]);
                if (found == secondary.boundary.nodes.end())
                    throw std::logic_error("HEX8 surface-contact slip output node mapping failed");
                const std::size_t output_node = static_cast<std::size_t>(found - secondary.boundary.nodes.begin());
                recovered_slip_weight[output_node] += contact_face_area[face_index];
                for (std::size_t component = 0; component < 3; ++component) {
                    recovered_total_slip[output_node][component] +=
                        contact_face_area[face_index] * contact_point_total_slip[face_index][local_node][component];
                    recovered_elastic_slip[output_node][component] +=
                        contact_face_area[face_index] * contact_point_elastic_slip[face_index][local_node][component];
                }
            }
        }
    }
    for (CartesianContactNodeSummary& summary : result)
        if (summary.tributary_area > 0.0) {
            summary.pressure = summary.contact_force / summary.tributary_area;
            summary.constraint_pressure = summary.pressure;
            summary.tangential_traction = summary.tangential_force / summary.tributary_area;
            const std::size_t node = static_cast<std::size_t>(&summary - result.data());
            if (!separately_averaged_friction && recovered_slip_weight[node] > 0.0)
                for (std::size_t component = 0; component < 3; ++component) {
                    summary.tangential_slip[component] =
                        recovered_total_slip[node][component] / recovered_slip_weight[node];
                    summary.elastic_tangential_slip[component] =
                        recovered_elastic_slip[node][component] / recovered_slip_weight[node];
                }
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
            if (node.constraint_pressure > 0.0) ++result.active_contact_nodes;
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
