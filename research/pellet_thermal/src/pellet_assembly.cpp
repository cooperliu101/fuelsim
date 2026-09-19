#include "core/thermal_assembly.hpp"
#include "dc3d8.hpp"
#include <algorithm>
#include <map>
#include <stdexcept>

namespace fuelsim::thermal {
std::vector<std::size_t> SpatialAssembly::pellet_surface_nodes(std::size_t r,
    const UnstructuredMeshMetadata& mesh) const {
    constexpr std::array<std::array<std::size_t, 4>, 6> faces{
        {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {0, 4, 7, 3}, {0, 3, 2, 1}, {4, 5, 6, 7}}};
    std::map<std::array<std::size_t, 4>, std::size_t> counts;
    for (std::size_t e = 0; e < _connectivity.size(); ++e)
        if (mesh.element_block_ids()[e] == _block_ids[r])
            for (const auto& face : faces) {
                std::array<std::size_t, 4> nodes{};
                for (std::size_t i = 0; i < 4; ++i)
                    nodes[i] = _connectivity[e][face[i]];
                std::sort(nodes.begin(), nodes.end());
                ++counts[nodes];
            }
    std::vector<std::size_t> boundary;
    for (const auto& entry : counts) {
        if (entry.second > 2)
            throw std::invalid_argument("Nonmanifold pellet mesh");
        if (entry.second == 1)
            boundary.insert(boundary.end(), entry.first.begin(), entry.first.end());
    }
    std::sort(boundary.begin(), boundary.end());
    boundary.erase(std::unique(boundary.begin(), boundary.end()), boundary.end());
    return boundary;
}

void SpatialAssembly::condense_pellets(const std::vector<std::size_t>& source_to_global) {
    if (std::all_of(_definition.regions.begin(), _definition.regions.end(), [](const RegionDefinition& region) {
            return region.pellet_response == "full";
        }))
        return;
    std::vector<Contribution> contributions;
    std::vector<std::size_t> sources;
    for (std::size_t r = 0; r < region_count(); ++r) {
        if (region(r).pellet_response == "full") {
            for (std::size_t e = 0; e < _contributions.size(); ++e)
                if (_contributions[e].region == r) {
                    contributions.push_back(std::move(_contributions[e]));
                    sources.push_back(_source_elements[e]);
                }
            continue;
        }
        const auto& material = region(r).material.functions->thermal;
        if (material.name != "constant_thermophysical")
            throw std::invalid_argument("Exact pellet condensation requires constant_thermophysical");
        std::vector<std::size_t> nodes;
        std::vector<std::vector<std::size_t>> region_connectivity;
        for (std::size_t e = 0; e < _contributions.size(); ++e)
            if (_contributions[e].region == r) {
                const auto& connectivity = _connectivity[_source_elements[e]];
                nodes.insert(nodes.end(), connectivity.begin(), connectivity.end());
                region_connectivity.push_back(connectivity);
            }
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        if (region(r).pellet_response == "surrogate") {
            Contribution pellet;
            pellet.region = r;
            pellet.pellet = _surrogates.size();
            _surrogates.emplace_back(region(r).pellet_model);
            _surrogates.back().validate_mesh(_source_nodes[r],
                _coordinates,
                elements::pellet_mesh_signature(_coordinates,
                    region_connectivity,
                    material.parameters.value("conductivity")));
            for (auto id : _source_nodes[r])
                pellet.dofs.push_back(source_to_global.at(id));
            contributions.push_back(std::move(pellet));
            sources.push_back(std::numeric_limits<std::size_t>::max());
            continue;
        }
        const auto n = nodes.size();
        std::vector<double> stiffness(n * n, 0.0), source(n, 0.0);
        const auto local_node = [&](std::size_t id) {
            return static_cast<std::size_t>(std::lower_bound(nodes.begin(), nodes.end(), id) - nodes.begin());
        };
        for (std::size_t e = 0; e < _contributions.size(); ++e) {
            const auto& c = _contributions[e];
            if (c.region != r)
                continue;
            const std::vector<double> temperature(8, region(r).initial_temperature), previous;
            const elements::ThermalInput
                input{material, c.geometry, temperature, previous, region(r).initial_temperature, 0.0, 0.0, 1.0};
            const auto result = elements::evaluate_dc3d8(input, true);
            const auto& connectivity = _connectivity[_source_elements[e]];
            for (std::size_t i = 0; i < 8; ++i) {
                const auto row = local_node(connectivity[i]);
                // Integrate unit volumetric source directly, avoiding cancellation
                // against a uniform-temperature conductivity residual.
                for (const auto& point : c.geometry.points)
                    source[row] += point.shape[i] * point.measure;
                for (std::size_t j = 0; j < 8; ++j)
                    stiffness[row * n + local_node(connectivity[j])] += result.jacobian[i * 8 + j];
            }
        }
        std::vector<std::size_t> boundary;
        Contribution pellet;
        pellet.region = r;
        pellet.pellet = _pellets.size();
        for (auto id : _source_nodes[r]) {
            boundary.push_back(local_node(id));
            pellet.dofs.push_back(source_to_global.at(id));
        }
        _pellets.emplace_back(stiffness, source, boundary);
        contributions.push_back(std::move(pellet));
        sources.push_back(std::numeric_limits<std::size_t>::max());
    }
    _contributions = std::move(contributions);
    _source_elements = std::move(sources);
}

elements::ThermalResult SpatialAssembly::evaluate_pellet(const Contribution& c,
    const std::vector<double>& temperature,
    double step,
    bool jacobian) const {
    if (step != 0.0)
        throw std::invalid_argument("Pellet responses support steady conduction only");
    elements::ThermalResult result;
    const double source = region_heat_source(c.region);
    if (region(c.region).pellet_response == "surrogate") {
        const auto& model = _surrogates.at(c.pellet);
        if (jacobian)
            model.evaluate_with_jacobian(temperature, source, result.residual, result.jacobian);
        else
            model.evaluate(temperature, source, result.residual);
        result.generated_heat_rate = source * model.volume();
    } else {
        const auto& pellet = _pellets.at(c.pellet);
        if (jacobian)
            pellet.evaluate_with_jacobian(temperature, source, result.residual, result.jacobian);
        else
            pellet.evaluate(temperature, source, result.residual);
        for (double value : pellet.source())
            result.generated_heat_rate += source * value;
    }
    return result;
}
} // namespace fuelsim::thermal
