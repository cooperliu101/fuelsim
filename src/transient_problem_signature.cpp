#include "transient_problem_signature.hpp"

#include "fuelsim/transient_problem.hpp"

#include <cstddef>
#include <cstring>
#include <string>

namespace fuelsim {
namespace {

constexpr std::uint64_t fnv_offset = 14695981039346656037ULL;
constexpr std::uint64_t fnv_prime = 1099511628211ULL;

void hash_bytes(std::uint64_t& hash, const void* data, std::size_t size) {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= static_cast<std::uint64_t>(bytes[index]);
        hash *= fnv_prime;
    }
}

void hash_size(std::uint64_t& hash, std::size_t value) {
    const std::uint64_t encoded = static_cast<std::uint64_t>(value);
    hash_bytes(hash, &encoded, sizeof(encoded));
}

void hash_integer(std::uint64_t& hash, std::int64_t value) {
    hash_bytes(hash, &value, sizeof(value));
}

void hash_double(std::uint64_t& hash, double value) {
    std::uint64_t encoded = 0;
    static_assert(sizeof(encoded) == sizeof(value));
    std::memcpy(&encoded, &value, sizeof(value));
    hash_bytes(hash, &encoded, sizeof(encoded));
}

void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_size(hash, value.size());
    hash_bytes(hash, value.data(), value.size());
}

void hash_thermoelastic(std::uint64_t& hash,
                        const ThermoelasticProperties& material) {
    hash_double(hash, material.conductivity_inverse_temperature);
    hash_double(hash, material.conductivity_offset);
    hash_double(hash, material.young_modulus);
    hash_double(hash, material.poisson_ratio);
    hash_double(hash, material.thermal_expansion);
    hash_double(hash, material.reference_temperature);
    hash_double(hash, material.young_modulus_temperature_coefficient);
    hash_double(hash, material.poisson_ratio_temperature_coefficient);
    hash_double(hash, material.thermal_expansion_temperature_coefficient);
}

} // namespace

std::uint64_t transient_problem_signature(const TransientProblem& problem) {
    std::uint64_t hash = fnv_offset;
    hash_size(hash, problem.dof_count());
    hash_size(hash, problem.region_count());
    const TransientProblemDefinition& definition = problem.definition();
    for (std::size_t region_value = 0; region_value < problem.region_count();
         ++region_value) {
        const RegionDefinition& spatial =
            definition.spatial.regions[region_value];
        const TransientInelasticProperties& transient =
            definition.regions[region_value].material;
        hash_string(hash, spatial.name);
        hash_string(hash, spatial.block);
        hash_integer(hash, spatial.block_id);
        hash_thermoelastic(hash, spatial.material);
        hash_double(hash, spatial.volumetric_heat_source);
        hash_double(hash, spatial.initial_temperature);
        hash_string(hash, spatial.heat_source_function);
        hash_integer(hash,
                     static_cast<std::int64_t>(spatial.strain_formulation));
        hash_double(hash, transient.density);
        hash_double(hash, transient.specific_heat);
        hash_integer(hash, static_cast<std::int64_t>(transient.behavior));
        hash_double(hash, transient.creep.coefficient);
        hash_double(hash, transient.creep.reference_stress);
        hash_double(hash, transient.creep.stress_exponent);
        hash_double(hash,
                    transient.creep.coefficient_temperature_coefficient);
        hash_double(hash,
                    transient.creep.reference_stress_temperature_coefficient);
        hash_double(hash,
                    transient.creep.stress_exponent_temperature_coefficient);
        hash_double(hash, transient.plasticity.yield_stress);
        hash_double(hash, transient.plasticity.isotropic_hardening_modulus);
        hash_double(
            hash,
            transient.plasticity.yield_stress_temperature_coefficient);
        hash_double(hash,
                    transient.plasticity.hardening_temperature_coefficient);

        const RegionMesh& mesh = problem.region_mesh(region_value);
        hash_size(hash, mesh.nodes().size());
        for (const RzPoint& point : mesh.nodes()) {
            hash_double(hash, point.r);
            hash_double(hash, point.z);
        }
        hash_size(hash, mesh.elements().size());
        for (const Quad4Element& element : mesh.elements()) {
            for (const std::size_t node : element.nodes)
                hash_size(hash, node);
        }
        for (const std::size_t source : mesh.source_node_ids())
            hash_size(hash, source);
        for (const std::size_t source : mesh.source_element_ids())
            hash_size(hash, source);
    }
    for (const ContactDefinition& contact : definition.spatial.contacts) {
        hash_string(hash, contact.name);
        hash_string(hash, contact.primary);
        hash_string(hash, contact.secondary);
        hash_integer(hash, contact.thermal ? 1 : 0);
        hash_integer(hash, contact.mechanical ? 1 : 0);
        hash_double(hash, contact.gap_conductivity);
        hash_double(hash, contact.minimum_gap);
        hash_double(hash, contact.penalty);
        hash_double(hash, contact.friction_coefficient);
        hash_integer(hash, contact.automatic_penalty ? 1 : 0);
        hash_double(hash, contact.penalty_factor);
        hash_integer(hash,
                     static_cast<std::int64_t>(
                         contact.mechanical_formulation));
        hash_double(hash, contact.penetration_tolerance);
        hash_size(hash, contact.maximum_augmented_iterations);
    }
    for (const BoundaryConditionDefinition& boundary :
         definition.spatial.boundary_conditions) {
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
        hash_integer(hash, boundary.use_displaced_geometry ? 1 : 0);
    }
    for (const PiecewiseLinearTimeTable& table :
         definition.spatial.time_tables) {
        hash_string(hash, table.name());
        hash_size(hash, table.times().size());
        for (std::size_t entry = 0; entry < table.times().size(); ++entry) {
            hash_double(hash, table.times()[entry]);
            hash_double(hash, table.values()[entry]);
        }
    }
    hash_size(hash, problem.contribution_count());
    for (std::size_t contribution = 0;
         contribution < problem.contribution_count(); ++contribution) {
        const LocalDofs dofs = problem.contribution_dofs(contribution);
        for (const std::size_t dof : dofs)
            hash_size(hash, dof);
    }
    return hash;
}

} // namespace fuelsim
