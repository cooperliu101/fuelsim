#include "io/problem_signature.hpp"
#include "core/problem_backend_access.hpp"
#include "detail/fnv_hash.hpp"
#include "fuelsim/core/transient_problem.hpp"
#include <cstddef>
#include <cstdint>
#include <string>

namespace fuelsim {
namespace {
void hash_size(std::uint64_t& hash, std::size_t value) {
    const std::uint64_t encoded = static_cast<std::uint64_t>(value);
    detail::fnv1a_bytes(hash, &encoded, sizeof(encoded));
}

void hash_integer(std::uint64_t& hash, std::int64_t value) { detail::fnv1a_bytes(hash, &value, sizeof(value)); }

void hash_double(std::uint64_t& hash, double value) {
    const std::uint64_t encoded = detail::encode_double_bits(value);
    detail::fnv1a_bytes(hash, &encoded, sizeof(encoded));
}

void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_size(hash, value.size());
    detail::fnv1a_bytes(hash, value.data(), value.size());
}

void hash_thermoelastic(std::uint64_t& hash, const ThermoelasticProperties& material) {
    hash_double(hash, material.reference_young_modulus);
    const std::uint64_t signature = material.functions->signature();
    detail::fnv1a_bytes(hash, &signature, sizeof(signature));
}

void hash_region_definition(std::uint64_t& hash, const RegionDefinition& spatial) {
    hash_string(hash, spatial.name);
    hash_string(hash, spatial.block);
    hash_integer(hash, spatial.block_id);
    hash_thermoelastic(hash, spatial.material);
    hash_double(hash, spatial.volumetric_heat_source);
    hash_double(hash, spatial.initial_temperature);
    hash_string(hash, spatial.heat_source_function);
    hash_integer(hash, static_cast<std::int64_t>(spatial.hex8_element_formulation));
}

void hash_boundaries(std::uint64_t& hash, const SpatialDefinition& definition, bool include_displaced_geometry) {
    for (const BoundaryConditionDefinition& boundary : definition.boundary_conditions) {
        hash_string(hash, boundary.name);
        hash_integer(hash, static_cast<std::int64_t>(boundary.type));
        hash_string(hash, boundary.boundary);
        hash_integer(hash, static_cast<std::int64_t>(boundary.field));
        hash_double(hash, boundary.value);
        hash_integer(hash, boundary.scale_with_load ? 1 : 0);
        hash_string(hash, boundary.function);
        hash_double(hash, boundary.heat_transfer_coefficient);
        hash_double(hash, boundary.ambient_temperature);
        hash_string(hash, boundary.coefficient_function);
        hash_string(hash, boundary.ambient_temperature_function);
        if (include_displaced_geometry) hash_integer(hash, boundary.use_displaced_geometry ? 1 : 0);
    }
}

void hash_time_tables(std::uint64_t& hash, const SpatialDefinition& definition) {
    for (const PiecewiseLinearTimeTable& table : definition.time_tables) {
        hash_string(hash, table.name());
        hash_size(hash, table.times().size());
        for (std::size_t entry = 0; entry < table.times().size(); ++entry) {
            hash_double(hash, table.times()[entry]);
            hash_double(hash, table.values()[entry]);
        }
    }
}

void hash_contacts(std::uint64_t& hash, const SpatialDefinition& definition) {
    for (const ContactDefinition& contact : definition.contacts) {
        hash_string(hash, contact.name);
        hash_string(hash, contact.primary);
        hash_string(hash, contact.secondary);
        hash_integer(hash, contact.thermal ? 1 : 0);
        hash_integer(hash, contact.mechanical ? 1 : 0);
        hash_double(hash, contact.gap_conductivity);
        hash_double(hash, contact.minimum_gap);
        hash_double(hash, contact.penalty);
        hash_double(hash, contact.friction_coefficient);
        hash_double(hash, contact.friction_slip_tolerance);
        hash_integer(hash, contact.automatic_penalty ? 1 : 0);
        hash_double(hash, contact.penalty_factor);
        hash_integer(hash, static_cast<std::int64_t>(contact.mechanical_formulation));
        hash_integer(hash, static_cast<std::int64_t>(contact.mechanical_discretization));
        hash_integer(hash, static_cast<std::int64_t>(contact.mechanical_sliding));
        hash_integer(hash, static_cast<std::int64_t>(contact.quad8_nodal_area_rule));
        hash_double(hash, contact.penetration_tolerance);
        hash_size(hash, contact.maximum_augmented_iterations);
    }
}
} // namespace

std::uint64_t transient_problem_signature(const TransientProblem& problem) {
    std::uint64_t hash = detail::fnv1a_offset;
    const bool cartesian = problem.is_cartesian_3d();
    const bool hex20 = cartesian && BackendAccess::cartesian_spatial(problem).uses_hex20();
    hash_string(hash, cartesian ? (hex20 ? "cartesian_3d_hex20_u2_t1" : "cartesian_3d_hex8") : "axisymmetric_rz_quad4");
    hash_string(hash, cartesian ? "xx,yy,zz,xy,yz,xz" : "rr,zz,hoop,rz");
    if (cartesian && hex20)
        hash_string(hash, "thermal_8_mechanical_27");
    else
        hash_size(hash, cartesian ? 8 : 4);
    hash_size(hash, problem.dof_count());
    if (cartesian) {
        const cartesian::SpatialAssembly& assembly = BackendAccess::cartesian_spatial(problem);
        const SpatialDefinition& definition = problem.definition();
        for (std::size_t region = 0; region < assembly.region_count(); ++region) {
            const RegionDefinition& spatial = definition.regions[region];
            hash_region_definition(hash, spatial);
            hash_integer(hash, static_cast<std::int64_t>(spatial.strain_formulation));
            if (hex20) {
                const Hex20RegionMesh& mesh = assembly.hex20_region_mesh(region);
                hash_size(hash, mesh.nodes().size());
                for (const CartesianPoint3& point : mesh.nodes()) {
                    hash_double(hash, point.x);
                    hash_double(hash, point.y);
                    hash_double(hash, point.z);
                }
                hash_size(hash, mesh.elements().size());
                for (const Hex20Element& element : mesh.elements())
                    for (const std::size_t node : element.nodes) hash_size(hash, node);
                continue;
            }
            const Hex8RegionMesh& mesh = assembly.region_mesh(region);
            hash_size(hash, mesh.nodes().size());
            for (const CartesianPoint3& point : mesh.nodes()) {
                hash_double(hash, point.x);
                hash_double(hash, point.y);
                hash_double(hash, point.z);
            }
            hash_size(hash, mesh.elements().size());
            for (const Hex8Element& element : mesh.elements())
                for (const std::size_t node : element.nodes) hash_size(hash, node);
        }
        hash_contacts(hash, definition);
        hash_boundaries(hash, definition, false);
        hash_time_tables(hash, definition);
        return hash;
    }
    const rz::TransientBackendView backend = BackendAccess::transient(problem);
    const SpatialDefinition& definition = problem.definition();
    for (std::size_t region = 0; region < backend.spatial.region_count(); ++region) {
        const RegionDefinition& spatial = definition.regions[region];
        hash_region_definition(hash, spatial);
        hash_integer(hash, static_cast<std::int64_t>(spatial.strain_formulation));
        const RegionMesh& mesh = backend.spatial.region_mesh(region);
        hash_size(hash, mesh.nodes().size());
        for (const RzPoint& point : mesh.nodes()) {
            hash_double(hash, point.r);
            hash_double(hash, point.z);
        }
        hash_size(hash, mesh.elements().size());
        for (const Quad4Element& element : mesh.elements())
            for (const std::size_t node : element.nodes) hash_size(hash, node);
    }
    hash_contacts(hash, definition);
    hash_boundaries(hash, definition, true);
    hash_time_tables(hash, definition);
    return hash;
}
} // namespace fuelsim
