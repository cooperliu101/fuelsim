#ifndef FUELSIM_TRANSIENT_FUEL_CLADDING_PROBLEM_HPP
#define FUELSIM_TRANSIENT_FUEL_CLADDING_PROBLEM_HPP

#include <cstddef>
#include <vector>

#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz_transient.hpp"
#include "fuelsim/steady_fuel_cladding_problem.hpp"

namespace fuelsim {

struct TransientFuelCladdingParameters final {
    SteadyFuelCladdingParameters steady;
    TransientInelasticProperties fuel;
    TransientInelasticProperties cladding;
};

struct TransientStepInput final {
    double end_time;
    double volumetric_heat_source;
};

struct RegionInelasticSummary final {
    double maximum_equivalent_plastic_strain;
    double maximum_equivalent_creep_strain;
};

class TransientFuelCladdingProblem final : public NonlinearProblem {
  public:
    explicit TransientFuelCladdingProblem(
        TransientFuelCladdingParameters parameters);
    TransientFuelCladdingProblem(TransientFuelCladdingParameters parameters,
                                 StructuredRzMesh fuel_mesh,
                                 StructuredRzMesh cladding_mesh);

    const TransientFuelCladdingParameters& parameters() const noexcept;
    const SteadyFuelCladdingProblem& steady_problem() const noexcept;
    const Quad4RzTransientKernel& fuel_kernel() const noexcept;
    const Quad4RzTransientKernel& cladding_kernel() const noexcept;

    const std::vector<double>& committed_solution() const noexcept;
    double committed_time() const noexcept;
    double committed_heat_source() const noexcept;
    bool time_step_active() const noexcept;
    double active_time_step() const;
    double active_end_time() const;

    void begin_time_step(const TransientStepInput& input);
    void commit_time_step(const std::vector<double>& converged_solution);
    void rollback_time_step() noexcept;

    const Quad4MaterialHistory&
    fuel_material_history(std::size_t element_index) const;
    const Quad4MaterialHistory&
    cladding_material_history(std::size_t element_index) const;
    RegionInelasticSummary summarize_fuel_history() const noexcept;
    RegionInelasticSummary summarize_cladding_history() const noexcept;

    InterfaceSummary
    summarize_interface(const std::vector<double>& state) const;

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

    TransientFuelCladdingParameters _parameters;
    SteadyFuelCladdingProblem _steady_problem;
    Quad4RzTransientKernel _fuel_kernel;
    Quad4RzTransientKernel _cladding_kernel;
    std::vector<Quad4MaterialHistory> _fuel_material_history;
    std::vector<Quad4MaterialHistory> _cladding_material_history;
    std::vector<double> _committed_solution;
    double _committed_time;
    double _committed_heat_source;
    double _active_time_step;
    double _active_end_time;
    double _active_heat_source;
    bool _time_step_active;
};

} // namespace fuelsim

#endif
