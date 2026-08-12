#ifndef FUELSIM_CARTESIAN3D_PROBLEM_ACCESS_HPP
#define FUELSIM_CARTESIAN3D_PROBLEM_ACCESS_HPP

#include "fuelsim/hex8.hpp"
#include "fuelsim/hex8_dof_map.hpp"
#include "fuelsim/hex8_thermoelastic.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim {
namespace cartesian3d {

struct TransientCommittedState final {
    std::vector<double> solution;
    std::vector<std::vector<std::array<SymmetricTensor3Values, 8>>> stresses;
    TransientConservationSummary conservation;
    double time = 0.0;
    double load_factor = 0.0;
};

class ProblemAccess final {
  public:
    static const Hex8DofMap& dof_map(const SteadyProblem& problem) noexcept;
    static std::size_t region_count(const SteadyProblem& problem) noexcept;
    static const RegionDefinition& region(const SteadyProblem& problem, std::size_t region_index);
    static const Hex8RegionMesh& region_mesh(const SteadyProblem& problem, std::size_t region_index);
    static const Hex8ThermoelasticKernel& region_kernel(const SteadyProblem& problem, std::size_t region_index);
    static std::size_t region_node_offset(const SteadyProblem& problem, std::size_t region_index);
    static const Hex8Geometry& region_element_geometry(const SteadyProblem& problem, std::size_t region_index, std::size_t element_index);

    static const TransientProblemDefinition& definition(const TransientProblem& problem) noexcept;
    static const Hex8DofMap& dof_map(const TransientProblem& problem) noexcept;
    static std::size_t region_count(const TransientProblem& problem) noexcept;
    static const RegionDefinition& region(const TransientProblem& problem, std::size_t region_index);
    static const Hex8RegionMesh& region_mesh(const TransientProblem& problem, std::size_t region_index);
    static const Hex8ThermoelasticKernel& region_kernel(const TransientProblem& problem, std::size_t region_index);
    static std::size_t region_node_offset(const TransientProblem& problem, std::size_t region_index);
    static const Hex8Geometry& region_element_geometry(const TransientProblem& problem, std::size_t region_index, std::size_t element_index);
    static const std::array<SymmetricTensor3Values, 8>& stress(const TransientProblem& problem, std::size_t region_index, std::size_t element_index);
    static TransientCommittedState committed_state(const TransientProblem& problem);
    static void restore_committed_state(TransientProblem& problem, TransientCommittedState state);
};

} // namespace cartesian3d
} // namespace fuelsim

#endif
