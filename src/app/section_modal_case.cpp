#include "app/section_modal_case.hpp"
#include "io/results_io.hpp"
#include "solver/section_modal_solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace fuelsim {
void run_section_modal_case(const FuelSimCaseDefinition& definition, const PetscSession& session) {
    if (definition.problem != CaseProblem::steady || definition.geometry != CaseGeometry::cartesian_3d
        || !definition.outputs.exodus_file.empty() || !definition.outputs.checkpoint_file.empty()
        || !definition.restart_file.empty() || definition.outputs.history_file.empty()
        || definition.outputs.csv_file.empty())
        throw std::invalid_argument("SectionModes requires steady cartesian_3d, csv and history outputs, without "
                                    "restart, checkpoint or Exodus results");
    const auto mesh = read_exodus_hex20(definition.mesh_file);
    if (definition.steady_execution.load_steps != 1)
        throw std::invalid_argument("SectionModes solves the linear final load directly and requires load_steps=1");
    const auto result = solve_section_modal(mesh,
        definition.spatial,
        definition.section_modes,
        definition.section_end_regions,
        definition.section_width_lines,
        definition.section_solid_ends,
        definition.solver);
    session.collective_root_action([&]() {
        std::ofstream summary(definition.outputs.csv_file), history(definition.outputs.history_file);
        if (!summary || !history)
            throw std::runtime_error("Cannot open section modal result files");
        summary << std::setprecision(17) << "metric,value\nmode_count," << result.mode_count << "\naxial_nodes,"
                << result.axial_nodes << "\npoisson_relaxation_modes," << result.poisson_modes
                << "\nautomatic_distortion_modes," << result.distortion_modes << "\nindependent_shear_modes,"
                << result.shear_modes << "\nindependent_axial_warping_modes," << result.axial_warping_modes
                << "\nshear_free_distortion_modes," << result.shear_free_modes << "\ntransverse_corrector_modes,"
                << result.transverse_corrector_modes << "\nrefinement_iterations," << result.refinement_iterations
                << "\nalgebraic_relative_residual," << result.algebraic_relative_residual << "\ndof_count,"
                << result.global_dof_count << "\nsection_basis_size," << result.basis_mode_count
                << "\nrecovered_local_dofs," << result.condensed_dof_count << "\nactive_axial_dofs,"
                << result.amplitudes.size() - result.solid_dof_count << "\nconstraint_count," << result.constraint_count
                << "\nstrain_energy," << result.energy << "\nequilibrium_relative_residual,"
                << result.equilibrium_relative_residual << "\nconstraint_maximum_error,"
                << result.constraint_maximum_error << "\nexternal_work," << result.external_work
                << "\nglobal_constraint_count," << result.global_constraint_count << "\nglobal_system_size,"
                << result.global_dof_count + result.global_constraint_count << "\nlocal_system_size,"
                << result.condensed_dof_count + result.constraint_count - result.global_constraint_count
                << "\nsection_preprocessing_seconds," << result.section_preprocessing_seconds << "\nassembly_seconds,"
                << result.assembly_seconds << "\ncondensation_seconds," << result.condensation_seconds
                << "\nsolve_and_refinement_seconds," << result.solve_and_refinement_seconds
                << "\nfield_recovery_seconds," << result.field_recovery_seconds << "\nwidth_line_count,"
                << result.width_line_count << "\nseed_mode_count," << result.seed_mode_count
                << "\nselected_refinement_iteration," << result.selected_refinement_iteration
                << "\nlast_refinement_relative_residual," << result.last_refinement_relative_residual
                << "\nwidth_enrichment_modes," << result.width_modes << "\nsolid_displacement_dofs,"
                << result.solid_dof_count << "\nsolid_elements," << result.solid_element_count << "\nmodal_elements,"
                << result.modal_element_count << "\nsolid_lower_interface," << result.solid_lower_interface
                << "\nsolid_upper_interface," << result.solid_upper_interface << '\n';
        // One inspectable result stream: fixed columns, record kind determines valid fields.
        history << std::setprecision(17)
                << "kind,id,x,y,z,ux,uy,uz,weight,exx,eyy,ezz,exy,eyz,exz,sxx,syy,szz,sxy,syz,sxz,"
                   "axial_force,moment_x,moment_y,torque,energy_per_length,mode,q,q_z,q_zz\n";
        for (std::size_t i = 0; i < mesh.nodes().size(); ++i) {
            const auto& p = mesh.nodes()[i];
            const auto& u = result.displacement[i];
            history << "node," << i << ',' << p.x << ',' << p.y << ',' << p.z << ',' << u.x << ',' << u.y << ',' << u.z
                    << ",0,0,0,0,0,0,0,0,0,0,0,0,0,,,,,,,,,\n";
        }
        for (std::size_t i = 0; i < result.samples.size(); ++i) {
            const auto& p = result.samples[i];
            history << (p.solid ? "solid_point," : "point,") << i << ',' << p.position.x << ',' << p.position.y << ','
                    << p.position.z << ',' << p.displacement.x << ',' << p.displacement.y << ',' << p.displacement.z
                    << ',' << p.weight;
            for (double value : p.strain)
                history << ',' << value;
            for (double value : p.stress)
                history << ',' << value;
            history << ",,,,,,,,,\n";
        }
        for (std::size_t i = 0; i < result.resultants.size(); ++i) {
            const auto& r = result.resultants[i];
            history << "section," << i << ",0,0," << r.z;
            for (std::size_t field = 0; field < 16; ++field)
                history << ",0";
            history << ',' << r.axial_force << ',' << r.moment_x << ',' << r.moment_y << ',' << r.torque << ','
                    << r.energy_per_length << ",,,,\n";
        }
        for (std::size_t node = 0; node < result.axial_nodes; ++node)
            for (std::size_t mode = 0; mode < result.nodal_mode_counts[node]; ++mode) {
                history << (result.nodal_mode_counts[node] > result.mode_count ? "local_amplitude," : "amplitude,")
                        << node << ",0,0," << result.axial_coordinates[node];
                for (std::size_t field = 0; field < 16; ++field)
                    history << ",0";
                history << ",,,,,," << mode;
                for (std::size_t order = 0; order < 3; ++order)
                    history << ','
                            << result.amplitudes[result.amplitude_map[(3 * mode + order) * result.axial_nodes + node]];
                history << '\n';
            }
        if (!summary || !history)
            throw std::runtime_error("Section modal result write failed");
        if (definition.outputs.console)
            std::cout << std::setprecision(12) << "section_modes=" << result.mode_count
                      << "\nmodal_dofs=" << result.global_dof_count
                      << "\nrecovered_local_dofs=" << result.condensed_dof_count << "\nstrain_energy=" << result.energy
                      << "\nequilibrium_relative_residual=" << result.equilibrium_relative_residual << '\n';
    });
}
} // namespace fuelsim
