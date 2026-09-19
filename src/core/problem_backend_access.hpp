#pragma once
#include "cartesian3d_assembly.hpp"
#include "cax2t_gps.hpp"
#include "cax4_types.hpp"
#include "contact_types.hpp"
#include "core/element_region_data.hpp"
#include "core/steady_problem.hpp"
#include "core/transient_problem.hpp"
#include "material_types.hpp"
#include "plane_assembly.hpp"
#include "radial_assembly.hpp"
#include "rz8_assembly.hpp"
#include "rz_assembly.hpp"
#include "thermal_assembly.hpp"
#include <array>
#include <vector>

namespace fuelsim {
enum class CommitPreparationStage;

struct TransientCommittedState final {
    std::vector<double> solution, previous_solution;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<CartesianMaterialHistory>> cartesian_material_histories;
    std::vector<std::vector<ContactPointHistory>> contact_histories;
    std::vector<double> raw_residual, external_load_residual;
    TransientConservationSummary conservation;
    double time = 0.0, load_factor = 0.0, previous_time = 0.0;
    std::vector<std::vector<Quad8MaterialHistory>> quad8_material_histories;
    std::vector<std::vector<Cax2tGpsMaterialHistory>> radial_material_histories;
};

namespace rz {
struct SteadyBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<AxisymmetricRegionData>& kernel_data;
};

struct TransientBackendView final {
    const SpatialAssembly& spatial;
    const std::vector<AxisymmetricRegionData>& kernel_data;
    const std::vector<std::vector<Quad4MaterialHistory>>& histories;
    const std::vector<double>& committed_solution;
    double active_time_step;
    bool time_step_active, include_thermal_time_term;
};
} // namespace rz

class BackendAccess final {
  public:
    static void set_commit_test_hook(TransientProblem& problem, std::function<void(CommitPreparationStage)> hook);
    static bool uses_thermal(const SteadyProblem& problem);
    static bool uses_thermal(const TransientProblem& problem);
    static const thermal::SpatialAssembly& thermal_spatial(const SteadyProblem& problem);
    static const thermal::SpatialAssembly& thermal_spatial(const TransientProblem& problem);
    static const plane::SpatialAssembly& plane_spatial(const SteadyProblem& problem);
    static const plane::SpatialAssembly& plane_spatial(const TransientProblem& problem);
    static bool uses_radial_gps(const SteadyProblem& problem);
    static bool uses_radial_gps(const TransientProblem& problem);
    static const radial::SpatialAssembly& radial_spatial(const SteadyProblem& problem);
    static const radial::SpatialAssembly& radial_spatial(const TransientProblem& problem);
    static const std::vector<std::vector<Cax2tGpsMaterialHistory>>& radial_material_histories(
        const TransientProblem& problem);
    static const rz8::SpatialAssembly& quad8_spatial(const SteadyProblem& problem);
    static const rz8::SpatialAssembly& quad8_spatial(const TransientProblem& problem);
    static bool uses_quad8(const SteadyProblem& problem);
    static const std::vector<AxisymmetricRegionData>& quad8_kernel_data(const SteadyProblem& problem);
    static const std::vector<std::vector<Quad8MaterialHistory>>& quad8_material_histories(
        const TransientProblem& problem);
    static rz::SteadyBackendView steady(const SteadyProblem& problem);
    static rz::TransientBackendView transient(const TransientProblem& problem);
    static const cartesian::SpatialAssembly& cartesian_spatial(const SteadyProblem& problem);
    static const cartesian::SpatialAssembly& cartesian_spatial(const TransientProblem& problem);
    static const std::vector<std::vector<CartesianMaterialHistory>>& cartesian_material_histories(
        const TransientProblem& problem);
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static const std::vector<double>& committed_raw_residual(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};
} // namespace fuelsim
