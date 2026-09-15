#pragma once
#include "core/modal_beam.hpp"
#include "core/spatial_definition.hpp"
#include "solver/petsc_solver.hpp"

namespace fuelsim {
struct SectionModalSample final {
    std::size_t axial_element, section_point;
    CartesianPoint3 position;
    double weight;
    SectionStrain strain, stress;
    CartesianPoint3 displacement;
};

struct SectionModalResultant final {
    double z = 0.0, axial_force = 0.0, moment_x = 0.0, moment_y = 0.0, torque = 0.0, energy_per_length = 0.0;
};

struct SectionModalResult final {
    std::vector<double> amplitudes;
    std::vector<CartesianPoint3> displacement;
    std::vector<SectionModalSample> samples;
    std::vector<SectionModalResultant> resultants;
    std::vector<double> axial_coordinates;
    std::vector<std::size_t> nodal_mode_counts, amplitude_map;
    std::size_t mode_count = 0, axial_nodes = 0, constraint_count = 0;
    std::size_t basis_mode_count = 0, global_dof_count = 0, condensed_dof_count = 0;
    std::size_t global_constraint_count = 0;
    std::size_t poisson_modes = 0, distortion_modes = 0, shear_modes = 0, axial_warping_modes = 0;
    std::size_t shear_free_modes = 0, transverse_corrector_modes = 0, refinement_iterations = 0;
    double energy = 0.0, equilibrium_relative_residual = 0.0, constraint_maximum_error = 0.0;
    double external_work = 0.0;
    double algebraic_relative_residual = 0.0;
    double section_preprocessing_seconds = 0.0, assembly_seconds = 0.0, condensation_seconds = 0.0;
    double solve_and_refinement_seconds = 0.0, field_recovery_seconds = 0.0;
};

// Small-strain, fixed extruded section. Original physical displacement and surface
// traction boundary conditions are projected; no section nodes are global unknowns.
SectionModalResult solve_section_modal(const UnstructuredHex20Mesh& source,
    const SpatialDefinition& definition,
    std::size_t mode_count,
    const std::vector<ModalEndRegion>& end_regions,
    const SolverOptions& options);
} // namespace fuelsim
