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
    const auto result = solve_section_modal(mesh, definition.spatial, definition.section_modes, definition.solver);
    session.collective_root_action([&]() {
        std::ofstream summary(definition.outputs.csv_file), history(definition.outputs.history_file);
        if (!summary || !history)
            throw std::runtime_error("Cannot open section modal result files");
        summary << std::setprecision(17) << "metric,value\nmode_count," << result.mode_count << "\naxial_nodes,"
                << result.axial_nodes << "\npoisson_relaxation_modes," << result.poisson_modes
                << "\nautomatic_distortion_modes," << result.distortion_modes << "\nindependent_shear_modes,"
                << result.shear_modes << "\ndof_count," << result.amplitudes.size() << "\nconstraint_count,"
                << result.constraint_count << "\nstrain_energy," << result.energy << "\nequilibrium_relative_residual,"
                << result.equilibrium_relative_residual << "\nconstraint_maximum_error,"
                << result.constraint_maximum_error << "\nexternal_work," << result.external_work << '\n';
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
            history << "point," << i << ',' << p.position.x << ',' << p.position.y << ',' << p.position.z << ",0,0,0,"
                    << p.weight;
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
            for (std::size_t mode = 0; mode < result.mode_count; ++mode) {
                history << "amplitude," << node << ",0,0," << result.axial_coordinates[node];
                for (std::size_t field = 0; field < 16; ++field)
                    history << ",0";
                history << ",,,,,," << mode;
                for (std::size_t order = 0; order < 3; ++order)
                    history << ',' << result.amplitudes[(3 * mode + order) * result.axial_nodes + node];
                history << '\n';
            }
        if (!summary || !history)
            throw std::runtime_error("Section modal result write failed");
        if (definition.outputs.console)
            std::cout << std::setprecision(12) << "section_modes=" << result.mode_count
                      << "\nmodal_dofs=" << result.amplitudes.size() << "\nstrain_energy=" << result.energy
                      << "\nequilibrium_relative_residual=" << result.equilibrium_relative_residual << '\n';
    });
}
} // namespace fuelsim
