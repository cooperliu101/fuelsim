#ifndef FUELSIM_TRANSIENT_PROBLEM_HPP
#define FUELSIM_TRANSIENT_PROBLEM_HPP

#include "fuelsim/inelastic_material.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz_transient.hpp"
#include "fuelsim/steady_problem.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fuelsim {

struct TransientRegionDefinition final {
    std::string region;
    TransientInelasticProperties material;
};

struct TransientProblemDefinition final {
    SteadyProblemDefinition spatial;
    std::vector<TransientRegionDefinition> regions;
};

struct TransientStepInput final {
    double end_time;
    double load_factor;
};

struct RegionInelasticSummary final {
    double maximum_equivalent_plastic_strain;
    double maximum_equivalent_creep_strain;
};

struct TransientCommittedState final {
    std::vector<double> solution;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>
        material_stresses;
    double time = 0.0;
    double load_factor = 0.0;
};

class TransientProblem final : public NonlinearProblem {
  public:
    TransientProblem(TransientProblemDefinition definition,
                     const UnstructuredQuad4Mesh& source_mesh);

    const TransientProblemDefinition& definition() const noexcept;
    const DofMap& dof_map() const noexcept;
    std::size_t region_count() const noexcept;
    std::size_t region_index(const std::string& name) const;
    std::size_t region_node_offset(std::size_t region_index) const;
    const RegionDefinition& region(std::size_t region_index) const;
    const RegionMesh& region_mesh(std::size_t region_index) const;
    const Quad4RzTransientKernel& region_kernel(std::size_t region_index) const;
    const Quad4RzGeometry&
    region_element_geometry(std::size_t region_index,
                            std::size_t element_index) const;

    const std::vector<double>& committed_solution() const noexcept;
    double committed_time() const noexcept;
    double committed_load_factor() const noexcept;
    bool time_step_active() const noexcept;
    double active_time_step() const;
    double active_end_time() const;
    std::vector<double> time_events() const;

    std::uint64_t committed_state_signature() const;
    TransientCommittedState committed_state() const;
    void restore_committed_state(TransientCommittedState state);

    void begin_time_step(const TransientStepInput& input);
    void commit_time_step(const std::vector<double>& converged_solution);
    void rollback_time_step() noexcept;

    const Quad4MaterialHistory&
    material_history(std::size_t region_index, std::size_t element_index) const;
    const std::array<AxisymmetricStressValues, 4>&
    material_stress(std::size_t region_index, std::size_t element_index) const;
    RegionInelasticSummary
    summarize_region_history(std::size_t region_index) const;
    InterfaceSummary
    summarize_interface(std::size_t contact_index,
                        const std::vector<double>& state) const;
    std::vector<ContactNodeSummary>
    summarize_contact_nodes(std::size_t contact_index,
                            const std::vector<double>& state) const;
    std::vector<std::size_t>
    contact_secondary_source_nodes(std::size_t contact_index) const;

    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    const std::vector<DirichletCondition>&
    dirichlet_conditions() const noexcept override;
    LocalDofs contribution_dofs(std::size_t contribution_index) const override;
    LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const LocalValues& state) const override;
    LocalSystem linearize_contribution(std::size_t contribution_index,
                                       const LocalValues& state) const override;

  protected:
    void add_state_independent_residual(
        std::vector<double>& residual) const override;

  private:
    Quad4TemperatureHistory
    committed_element_temperature(std::size_t contribution_index) const;
    void require_active_time_step() const;

    TransientProblemDefinition _definition;
    SteadyProblem _spatial_model;
    std::vector<Quad4RzTransientKernel> _region_kernels;
    std::vector<std::vector<Quad4MaterialHistory>> _material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>
        _material_stresses;
    std::vector<double> _committed_solution;
    double _committed_time;
    double _committed_load_factor;
    double _active_time_step;
    double _active_end_time;
    double _active_load_factor;
    bool _time_step_active;
};

} // namespace fuelsim

#endif
