#pragma once
#include "cartesian3d_assembly.hpp"
#include "cax4_types.hpp"
#include "contact_types.hpp"
#include "core/element_region_data.hpp"
#include "core/steady_problem.hpp"
#include "core/transient_problem.hpp"
#include "material_types.hpp"
#include "rz8_assembly.hpp"
#include "rz_assembly.hpp"
#include <array>
#include <vector>

namespace fuelsim {
struct TransientCommittedState final {
    std::vector<double> solution, previous_solution;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<CartesianMaterialHistory>> cartesian_material_histories;
    std::vector<std::vector<ContactPointHistory>> contact_histories;
    std::vector<double> raw_residual, external_load_residual;
    TransientConservationSummary conservation;
    double time = 0.0, load_factor = 0.0, previous_time = 0.0;
    std::vector<std::vector<Quad8MaterialHistory>> quad8_material_histories;
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
    static const rz8::SpatialAssembly& quad8_spatial(const SteadyProblem& problem) noexcept;
    static const rz8::SpatialAssembly& quad8_spatial(const TransientProblem& problem) noexcept;
    static bool uses_quad8(const SteadyProblem& problem) noexcept;
    static const std::vector<AxisymmetricRegionData>& quad8_kernel_data(const SteadyProblem& problem) noexcept;
    static const std::vector<std::vector<Quad8MaterialHistory>>& quad8_material_histories(
        const TransientProblem& problem) noexcept;
    static rz::SteadyBackendView steady(const SteadyProblem& problem) noexcept;
    static rz::TransientBackendView transient(const TransientProblem& problem) noexcept;
    static const cartesian::SpatialAssembly& cartesian_spatial(const SteadyProblem& problem) noexcept;
    static const cartesian::SpatialAssembly& cartesian_spatial(const TransientProblem& problem) noexcept;
    static const std::vector<std::vector<CartesianMaterialHistory>>& cartesian_material_histories(
        const TransientProblem& problem) noexcept;
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static const std::vector<double>& committed_raw_residual(const TransientProblem& problem) noexcept;
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};
} // namespace fuelsim
