#include "core/cross_section.hpp"
#include "quad8_face.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fuelsim {
CrossSection::CrossSection(std::vector<CartesianPoint3> nodes,
    std::vector<SectionCell> cells,
    std::vector<SectionRegion> regions)
    : _nodes(std::move(nodes)), _regions(std::move(regions)) {
    if (_nodes.empty() || cells.empty() || _regions.empty())
        throw std::invalid_argument("CrossSection requires nodes, cells and material regions");
    for (const auto& node : _nodes)
        if (!std::isfinite(node.x) || !std::isfinite(node.y) || node.z != 0.0)
            throw std::invalid_argument("CrossSection requires finite coordinates in the xy plane");
    for (std::size_t i = 0; i < _regions.size(); ++i) {
        const auto& region = _regions[i];
        const auto& functions = region.material.functions();
        if (region.name.empty() || !std::isfinite(region.temperature) || functions.has_creep()
            || functions.has_plasticity() || !functions.eigenstrains.empty())
            throw std::invalid_argument("CrossSection requires named, purely elastic regions without eigenstrain");
        for (std::size_t j = 0; j < i; ++j)
            if (region.name == _regions[j].name)
                throw std::invalid_argument("CrossSection region names must be unique");
    }
    std::vector<std::vector<std::size_t>> adjacency(_nodes.size());
    std::vector<bool> region_used(_regions.size(), false);
    double elastic_area = 0.0;
    for (const auto& cell : cells) {
        if (cell.region >= _regions.size())
            throw std::invalid_argument("CrossSection cell region is out of range");
        region_used[cell.region] = true;
        Quad8FaceCoordinates coordinates{};
        for (std::size_t i = 0; i < 8; ++i) {
            const auto id = cell.element.nodes[i];
            if (id >= _nodes.size())
                throw std::invalid_argument("CrossSection connectivity is out of range");
            for (std::size_t j = 0; j < i; ++j)
                if (id == cell.element.nodes[j])
                    throw std::invalid_argument("CrossSection cell has repeated nodes");
            coordinates[i] = _nodes[id];
            adjacency[id].insert(adjacency[id].end(), cell.element.nodes.begin(), cell.element.nodes.end());
        }
        // A surface norm loses orientation: explicitly check the signed planar map.
        for (double xi : {-1.0, 0.0, 1.0})
            for (double eta : {-1.0, 0.0, 1.0}) {
                const auto g = make_quad8_face_mechanical_point(coordinates, xi, eta, 1.0);
                const double det = g.tangent_xi.x * g.tangent_eta.y - g.tangent_xi.y * g.tangent_eta.x;
                if (!std::isfinite(det) || det <= 0.0)
                    throw std::invalid_argument("CrossSection requires a positive planar Jacobian");
            }
        for (const auto& g : make_quad8_face_geometry(coordinates).mechanical_points) {
            const double det = g.tangent_xi.x * g.tangent_eta.y - g.tangent_xi.y * g.tangent_eta.x;
            if (!std::isfinite(det) || det <= 0.0)
                throw std::invalid_argument("CrossSection quadrature Jacobian is invalid");
            SectionPoint point{};
            point.nodes = cell.element.nodes;
            point.shape = g.displacement_shape;
            point.weight = det * g.quadrature_weight;
            point.region = cell.region;
            for (std::size_t i = 0; i < 8; ++i) {
                point.gradient[i] = {(g.derivative_xi[i] * g.tangent_eta.y - g.derivative_eta[i] * g.tangent_xi.y)
                                         / det,
                    (-g.derivative_xi[i] * g.tangent_eta.x + g.derivative_eta[i] * g.tangent_xi.x) / det};
                point.position.x += point.shape[i] * coordinates[i].x;
                point.position.y += point.shape[i] * coordinates[i].y;
            }
            const auto& region = _regions[point.region];
            const auto properties =
                region.material.active_properties(region.temperature, material_context(0.0, point.position));
            const double lambda = properties.lame_lambda.value(), mu = properties.shear_modulus.value();
            const double young = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
            if (!std::isfinite(young) || young <= 0.0 || mu <= 0.0 || 3.0 * lambda + 2.0 * mu <= 0.0)
                throw std::invalid_argument("CrossSection requires stable isotropic elasticity");
            _area += point.weight;
            _centroid.x += point.weight * point.position.x;
            _centroid.y += point.weight * point.position.y;
            elastic_area += point.weight * young;
            _elastic_center.x += point.weight * young * point.position.x;
            _elastic_center.y += point.weight * young * point.position.y;
            _points.push_back(point);
        }
    }
    std::vector<bool> visited(_nodes.size(), false);
    std::vector<std::size_t> pending{0};
    visited[0] = true;
    for (std::size_t i = 0; i < pending.size(); ++i)
        for (auto neighbor : adjacency[pending[i]])
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                pending.push_back(neighbor);
            }
    if (pending.size() != _nodes.size()
        || std::find(region_used.begin(), region_used.end(), false) != region_used.end())
        throw std::invalid_argument("CrossSection requires one connected mesh and no unused nodes or regions");
    _centroid.x /= _area;
    _centroid.y /= _area;
    _elastic_center.x /= elastic_area;
    _elastic_center.y /= elastic_area;
}

SectionKinematics classic_section_kinematics(const CrossSection& section,
    const ClassicSectionBasis& basis,
    std::size_t index,
    const std::array<double, 3>& q,
    const std::array<double, 3>& first,
    const std::array<double, 3>& second,
    const std::array<double, 3>& third) {
    for (const auto* values : {&q, &first, &second, &third})
        for (double value : *values)
            if (!std::isfinite(value))
                throw std::invalid_argument("Classic section amplitudes must be finite");
    if (!std::isfinite(basis.origin.x) || !std::isfinite(basis.origin.y))
        throw std::invalid_argument("Classic section origin must be finite");
    const auto& point = section.points().at(index);
    const std::size_t n = section.nodes().size();
    const double x = point.position.x - basis.origin.x, y = point.position.y - basis.origin.y;
    SectionKinematics result;
    result.displacement = {q[1], q[2], q[0] - x * first[1] - y * first[2]};
    result.gradient[2] = first[1];
    result.gradient[5] = first[2];
    result.gradient[6] = -first[1];
    result.gradient[7] = -first[2];
    result.gradient[8] = first[0] - x * second[1] - y * second[2];
    const std::array<double, 3> amplitude{first[0], second[1], second[2]};
    const std::array<double, 3> derivative{second[0], third[1], third[2]};
    for (std::size_t mode = 0; mode < 3; ++mode) {
        const auto& correction = basis.modes[mode].correction;
        if (correction.size() != 2 * n || static_cast<std::size_t>(basis.modes[mode].kind) != mode)
            throw std::invalid_argument("Classic section correction size or mode order mismatch");
        for (std::size_t i = 0; i < 8; ++i) {
            const double ux = correction[point.nodes[i]], uy = correction[n + point.nodes[i]];
            if (!std::isfinite(ux) || !std::isfinite(uy))
                throw std::invalid_argument("Classic section corrections must be finite");
            result.displacement[0] += point.shape[i] * ux * amplitude[mode];
            result.displacement[1] += point.shape[i] * uy * amplitude[mode];
            result.gradient[0] += point.gradient[i][0] * ux * amplitude[mode];
            result.gradient[1] += point.gradient[i][1] * ux * amplitude[mode];
            result.gradient[3] += point.gradient[i][0] * uy * amplitude[mode];
            result.gradient[4] += point.gradient[i][1] * uy * amplitude[mode];
            result.gradient[2] += point.shape[i] * ux * derivative[mode];
            result.gradient[5] += point.shape[i] * uy * derivative[mode];
        }
    }
    return result;
}

SectionStrain classic_section_strain(const CrossSection& section,
    const ClassicSectionBasis& basis,
    std::size_t index,
    const std::array<double, 3>& amplitude,
    const std::array<double, 3>& derivative) {
    const auto kinematics = classic_section_kinematics(section,
        basis,
        index,
        {},
        {amplitude[0], 0.0, 0.0},
        {derivative[0], amplitude[1], amplitude[2]},
        {0.0, derivative[1], derivative[2]});
    const auto& g = kinematics.gradient;
    return {g[0], g[4], g[8], 0.5 * (g[1] + g[3]), 0.5 * (g[5] + g[7]), 0.5 * (g[2] + g[6])};
}

SectionResponse evaluate_classic_section(const CrossSection& section,
    const ClassicSectionBasis& basis,
    const std::array<double, 3>& amplitude) {
    SectionResponse result;
    constexpr std::array<double, 6> metric{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    for (std::size_t q = 0; q < section.points().size(); ++q) {
        const auto& point = section.points()[q];
        const auto& region = section.regions()[point.region];
        const auto strain = classic_section_strain(section, basis, q, amplitude);
        const auto response = evaluate_stress_tangent(region.material,
            strain,
            region.temperature,
            0.0,
            nullptr,
            material_context(0.0, point.position));
        const SectionStrain stress{response.stress.xx,
            response.stress.yy,
            response.stress.zz,
            response.stress.xy,
            response.stress.yz,
            response.stress.xz};
        result.strain.push_back(strain);
        result.stress.push_back(stress);
        std::array<SectionStrain, 3> b{};
        for (std::size_t m = 0; m < 3; ++m) {
            std::array<double, 3> unit{};
            unit[m] = 1.0;
            b[m] = classic_section_strain(section, basis, q, unit);
        }
        for (std::size_t i = 0; i < 6; ++i) {
            result.energy += 0.5 * point.weight * metric[i] * strain[i] * stress[i];
            for (std::size_t m = 0; m < 3; ++m) {
                result.force[m] += point.weight * metric[i] * b[m][i] * stress[i];
                for (std::size_t k = 0; k < 3; ++k)
                    for (std::size_t j = 0; j < 6; ++j)
                        result.tangent[3 * m + k] +=
                            point.weight * metric[i] * b[m][i] * response.tangent[i][j] * b[k][j];
            }
        }
    }
    return result;
}
} // namespace fuelsim
