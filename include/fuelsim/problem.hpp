#ifndef FUELSIM_PROBLEM_HPP
#define FUELSIM_PROBLEM_HPP

#include <array>
#include <cstddef>
#include <vector>

#include "fuelsim/dof_map.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz.hpp"

namespace fuelsim {

struct M0Parameters final {
    double inner_radius;
    double outer_radius;
    double length;
    std::size_t radial_elements;
    std::size_t axial_elements;

    ThermoelasticProperties fuel;
    double volumetric_heat_source;
    double outer_temperature;
    double initial_temperature;

    double inner_pressure;
    double outer_pressure;
};

class M0Problem final : public NonlinearProblem {
  public:
    explicit M0Problem(M0Parameters parameters);

    const M0Parameters& parameters() const noexcept;
    const StructuredRzMesh& mesh() const noexcept;
    const DofMap& dof_map() const noexcept;
    const Quad4RzThermoelasticKernel& kernel() const noexcept;

    std::size_t dof_count() const noexcept override;
    std::size_t element_count() const noexcept;
    std::size_t contribution_count() const noexcept override;

    const std::vector<DirichletCondition>&
    dirichlet_conditions() const noexcept override;

    std::vector<double> initial_state() const;

    LocalDofs element_dofs(std::size_t element_index) const;
    LocalValues element_state(std::size_t element_index,
                              const std::vector<double>& global_state) const;
    LocalResidual element_residual(std::size_t element_index,
                                   const LocalValues& state) const;
    LocalSystem linearize_element(std::size_t element_index,
                                  const LocalValues& state) const;

    LocalDofs contribution_dofs(std::size_t contribution_index) const override;
    LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const LocalValues& state) const override;
    LocalSystem linearize_contribution(std::size_t contribution_index,
                                       const LocalValues& state) const override;

    const Quad4RzGeometry& element_geometry(std::size_t element_index) const;

  protected:
    void add_state_independent_residual(
        std::vector<double>& residual) const override;

  private:
    void build_dirichlet_conditions();
    void add_pressure_residual(std::vector<double>& residual) const;

    M0Parameters _parameters;
    StructuredRzMesh _mesh;
    DofMap _dof_map;
    Quad4RzThermoelasticKernel _kernel;
    std::vector<Quad4RzGeometry> _geometries;
    std::vector<DirichletCondition> _dirichlet_conditions;
};

} // namespace fuelsim

#endif
