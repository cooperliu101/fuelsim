#include "io/problem_signature.hpp"
#include "core/problem_backend_access.hpp"
#include "core/transient_problem.hpp"
#include "fnv_hash.hpp"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fuelsim {
namespace {
void hash_size(std::uint64_t& hash, std::size_t value) {
    const std::uint64_t encoded = static_cast<std::uint64_t>(value);
    hashing::fnv1a_bytes(hash, &encoded, sizeof(encoded));
}

void hash_integer(std::uint64_t& hash, std::int64_t value) {
    hashing::fnv1a_bytes(hash, &value, sizeof(value));
}

void hash_double(std::uint64_t& hash, double value) {
    const std::uint64_t encoded = hashing::encode_double_bits(value);
    hashing::fnv1a_bytes(hash, &encoded, sizeof(encoded));
}

void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_size(hash, value.size());
    hashing::fnv1a_bytes(hash, value.data(), value.size());
}

void hash_cartesian_nodes(std::uint64_t& hash, const std::vector<CartesianPoint3>& nodes) {
    hash_size(hash, nodes.size());
    for (const CartesianPoint3& point : nodes) {
        hash_double(hash, point.x);
        hash_double(hash, point.y);
        hash_double(hash, point.z);
    }
}

void hash_thermoelastic(std::uint64_t& hash, const ThermoelasticProperties& material) {
    hash_double(hash, material.reference_young_modulus);
    const std::uint64_t signature = material.functions->signature();
    hashing::fnv1a_bytes(hash, &signature, sizeof(signature));
}

void hash_region_definition(std::uint64_t& hash, const RegionDefinition& spatial) {
    hash_string(hash, spatial.name);
    hash_string(hash, spatial.block);
    hash_integer(hash, spatial.block_id);
    hash_thermoelastic(hash, spatial.material);
    hash_double(hash, spatial.volumetric_heat_source);
    hash_double(hash, spatial.initial_temperature);
    for (double acceleration : spatial.body_acceleration)
        hash_double(hash, acceleration);
    hash_string(hash, spatial.heat_source_function);
    hash_integer(hash, static_cast<std::int64_t>(spatial.heat_source_time_evaluation));
    hash_integer(hash, static_cast<std::int64_t>(spatial.hex8_element_formulation));
    hash_integer(hash, static_cast<std::int64_t>(spatial.hex20_element_formulation));
    hash_integer(hash, static_cast<std::int64_t>(spatial.rz_element_formulation));
    hash_integer(hash, static_cast<std::int64_t>(spatial.strain_formulation));
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
        if (include_displaced_geometry)
            hash_integer(hash, boundary.use_displaced_geometry ? 1 : 0);
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
        hash_integer(hash, static_cast<std::int64_t>(contact.thermal_discretization));
        hash_integer(hash, static_cast<std::int64_t>(contact.mechanical_sliding));
        hash_integer(hash, static_cast<std::int64_t>(contact.quad8_nodal_area_rule));
        hash_double(hash, contact.penetration_tolerance);
        hash_size(hash, contact.maximum_augmented_iterations);
    }
}
} // namespace

std::uint64_t transient_problem_signature(const TransientProblem& problem) {
    std::uint64_t hash = hashing::fnv1a_offset;
    if (BackendAccess::uses_thermal(problem)) {
        hash_string(hash, "thermal_fixed_order_dependent_capacity_v1");
        const auto& spatial = BackendAccess::thermal_spatial(problem);
        hash_integer(hash, spatial.axisymmetric() ? 1 : 0);
        for (const auto& r : problem.definition().regions) {
            hash_region_definition(hash, r);
            hash_integer(hash, static_cast<std::int64_t>(r.thermal_element));
        }
        for (const auto& point : spatial.coordinates()) {
            hash_double(hash, point.x);
            hash_double(hash, point.y);
            hash_double(hash, point.z);
        }
        for (std::size_t i = 0; i < spatial.contribution_count(); ++i) {
            std::vector<std::size_t> dofs;
            spatial.contribution_dofs(i, dofs);
            hash_size(hash, dofs.size());
            for (auto n : dofs)
                hash_size(hash, n);
        }
        hash_boundaries(hash, problem.definition(), true);
        hash_contacts(hash, problem.definition());
        for (const auto& contact : problem.definition().contacts) {
            hash_integer(hash, static_cast<std::int64_t>(contact.gap_heat_conductance_law));
            hash_double(hash, contact.gap_conductance);
            hash_double(hash, contact.gap_conductance_clearance_derivative);
            hash_double(hash, contact.gap_conductance_pressure_derivative);
            hash_double(hash, contact.gap_conductance_temperature_derivative);
            hash_double(hash, contact.gap_conductance_reference_temperature);
        }
        hash_time_tables(hash, problem.definition());
        return hash;
    }
    const bool cartesian = problem.is_cartesian_3d();
    if (problem.uses_plane_quad8()) {
        hash_string(hash, "generalized_plane_strain_cpeg8t_segmented_thermal_v3");
        hash_size(hash, problem.dof_count());
        const auto& spatial = BackendAccess::plane_spatial(problem);
        hash_size(hash, spatial.section_count());
        for (const auto& section : problem.definition().generalized_plane_strain) {
            hash_string(hash, section.name);
            hash_size(hash, section.blocks.size());
            for (const auto& block : section.blocks)
                hash_string(hash, block);
            hash_double(hash, section.initial_thickness);
            for (std::size_t i = 0; i < 3; ++i) {
                hash_integer(hash, section.prescribed[i] ? 1 : 0);
                if (section.prescribed[i])
                    hash_double(hash, *section.prescribed[i]);
                hash_string(hash, section.functions[i]);
            }
        }
        for (std::size_t r = 0; r < spatial.region_count(); ++r) {
            hash_region_definition(hash, spatial.region(r));
            hash_size(hash, spatial.region_element_count(r));
            for (std::size_t e = 0; e < spatial.region_element_count(r); ++e) {
                const auto& geometry = spatial.geometry(r, e);
                for (const auto& point : geometry.coordinates)
                    for (auto coordinate : point)
                        hash_double(hash, coordinate);
                hash_double(hash, geometry.thickness);
                for (auto coordinate : geometry.reference_point)
                    hash_double(hash, coordinate);
                std::vector<std::size_t> dofs;
                spatial.contribution_dofs(spatial.region_element_offset(r) + e, dofs);
                for (auto dof : dofs)
                    hash_size(hash, dof);
            }
        }
        for (const auto& condition : spatial.dirichlet_conditions())
            hash_size(hash, condition.dof);
        for (std::size_t index = spatial.volume_contribution_count(); index < spatial.contribution_count(); ++index) {
            std::vector<std::size_t> dofs;
            spatial.contribution_dofs(index, dofs);
            hash_size(hash, dofs.size());
            for (auto dof : dofs)
                hash_size(hash, dof);
            if (index < spatial.contact_offset()) {
                const auto identity = spatial.boundary_identity(index);
                hash_size(hash, identity.first);
                hash_size(hash, identity.second);
            }
        }
        hash_contacts(hash, problem.definition());
        for (const auto& contact : problem.definition().contacts) {
            hash_integer(hash, static_cast<std::int64_t>(contact.gap_heat_conductance_law));
            hash_double(hash, contact.gap_conductance);
            hash_double(hash, contact.gap_conductance_clearance_derivative);
            hash_double(hash, contact.gap_conductance_pressure_derivative);
            hash_double(hash, contact.gap_conductance_temperature_derivative);
            hash_double(hash, contact.gap_conductance_reference_temperature);
        }
        for (const auto& boundary : problem.definition().boundary_conditions)
            hash_integer(hash, boundary.configuration_explicit ? 1 : 0);
        hash_boundaries(hash, problem.definition(), true);
        hash_time_tables(hash, problem.definition());
        return hash;
    }
    if (BackendAccess::uses_radial_gps(problem)) {
        hash_string(hash, "axisymmetric_1d_bar2_gps_two_points");
        hash_size(hash, problem.dof_count());
        const auto& spatial = BackendAccess::radial_spatial(problem);
        const auto& mesh = spatial.source_mesh();
        hash_size(hash, mesh.nodes().size());
        for (const auto& node : mesh.nodes()) {
            hash_double(hash, node.r);
            hash_double(hash, node.z);
        }
        hash_size(hash, mesh.elements().size());
        for (std::size_t i = 0; i < mesh.elements().size(); ++i) {
            for (const auto node : mesh.elements()[i].nodes)
                hash_size(hash, node);
            for (const auto node : mesh.elements()[i].axial_nodes)
                hash_size(hash, node);
            hash_integer(hash, mesh.element_block_ids()[i]);
        }
        for (const auto& block : mesh.element_blocks()) {
            hash_integer(hash, block.id);
            hash_string(hash, block.name);
        }
        for (const auto& set : mesh.node_sets()) {
            hash_integer(hash, set.id);
            hash_string(hash, set.name);
            hash_size(hash, set.nodes.size());
            for (const auto node : set.nodes)
                hash_size(hash, node);
        }
        for (const auto& set : mesh.side_sets()) {
            hash_integer(hash, set.id);
            hash_string(hash, set.name);
            hash_size(hash, set.sides.size());
            for (const auto& side : set.sides) {
                hash_size(hash, side.element);
                hash_size(hash, side.local_side);
            }
        }
        for (const auto& region : problem.definition().regions)
            hash_region_definition(hash, region);
        hash_contacts(hash, problem.definition());
        for (const auto& contact : problem.definition().contacts) {
            hash_integer(hash, static_cast<std::int64_t>(contact.gap_heat_conductance_law));
            hash_double(hash, contact.gap_conductance);
            hash_double(hash, contact.gap_conductance_clearance_derivative);
            hash_double(hash, contact.gap_conductance_pressure_derivative);
            hash_double(hash, contact.gap_conductance_temperature_derivative);
            hash_double(hash, contact.gap_conductance_reference_temperature);
        }
        hash_boundaries(hash, problem.definition(), true);
        hash_time_tables(hash, problem.definition());
        return hash;
    }
    if (problem.uses_quad8()) {
        hash_string(hash, "axisymmetric_rz_quad8_u2_t1_nine_points");
        hash_size(hash, problem.dof_count());
        const auto& spatial = BackendAccess::quad8_spatial(problem);
        for (std::size_t r = 0; r < spatial.region_count(); ++r) {
            hash_region_definition(hash, spatial.region(r));
            const auto& mesh = spatial.region_mesh(r);
            hash_size(hash, mesh.nodes().size());
            for (std::size_t n = 0; n < mesh.nodes().size(); ++n) {
                hash_double(hash, mesh.nodes()[n].r);
                hash_double(hash, mesh.nodes()[n].z);
                hash_size(hash, mesh.source_node_ids()[n]);
                hash_size(hash, spatial.global_node(r, n));
                hash_integer(hash, mesh.temperature_nodes()[n] ? 1 : 0);
            }
            hash_size(hash, mesh.elements().size());
            for (const auto& e : mesh.elements())
                for (auto n : e.nodes)
                    hash_size(hash, n);
        }
        for (const auto& bc : problem.definition().boundary_conditions)
            hash_integer(hash, bc.configuration_explicit ? 1 : 0);
        for (const auto& bc : spatial.dirichlet_conditions())
            hash_size(hash, bc.dof);
        for (std::size_t i = spatial.volume_contribution_count(); i < spatial.sparsity_contribution_count(); ++i) {
            std::vector<std::size_t> dofs;
            spatial.sparsity_contribution_dofs(i, dofs);
            hash_size(hash, dofs.size());
            for (auto d : dofs)
                hash_size(hash, d);
        }
        hash_contacts(hash, problem.definition());
        hash_boundaries(hash, problem.definition(), true);
        hash_time_tables(hash, problem.definition());
        return hash;
    }
    const bool hex20 = cartesian && BackendAccess::cartesian_spatial(problem).uses_hex20();
    hash_string(hash, cartesian ? (hex20 ? "cartesian_3d_hex20_u2_t1" : "cartesian_3d_hex8") : "axisymmetric_rz_quad4");
    hash_string(hash, cartesian ? "xx,yy,zz,xy,yz,xz" : "rr,zz,hoop,rz");
    if (cartesian && hex20)
        hash_string(hash, "hex20_region_selected_quadrature");
    else
        hash_size(hash, cartesian ? 8 : 4);
    hash_size(hash, problem.dof_count());
    if (cartesian) {
        const cartesian::SpatialAssembly& assembly = BackendAccess::cartesian_spatial(problem);
        const SpatialDefinition& definition = problem.definition();
        for (std::size_t region = 0; region < assembly.region_count(); ++region) {
            const RegionDefinition& spatial = definition.regions[region];
            hash_region_definition(hash, spatial);
            if (hex20) {
                const Hex20RegionMesh& mesh = assembly.hex20_region_mesh(region);
                hash_cartesian_nodes(hash, mesh.nodes());
                hash_size(hash, mesh.elements().size());
                for (const Hex20Element& element : mesh.elements())
                    for (const std::size_t node : element.nodes)
                        hash_size(hash, node);
                continue;
            }
            const Hex8RegionMesh& mesh = assembly.region_mesh(region);
            hash_cartesian_nodes(hash, mesh.nodes());
            hash_size(hash, mesh.elements().size());
            for (const Hex8Element& element : mesh.elements())
                for (const std::size_t node : element.nodes)
                    hash_size(hash, node);
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
        const RegionMesh& mesh = backend.spatial.region_mesh(region);
        hash_size(hash, mesh.nodes().size());
        for (const RzPoint& point : mesh.nodes()) {
            hash_double(hash, point.r);
            hash_double(hash, point.z);
        }
        hash_size(hash, mesh.elements().size());
        for (const Quad4Element& element : mesh.elements())
            for (const std::size_t node : element.nodes)
                hash_size(hash, node);
    }
    hash_contacts(hash, definition);
    hash_boundaries(hash, definition, true);
    hash_time_tables(hash, definition);
    return hash;
}
} // namespace fuelsim
