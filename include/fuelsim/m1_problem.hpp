#ifndef FUELSIM_M1_PROBLEM_HPP
#define FUELSIM_M1_PROBLEM_HPP

#include <cstddef>
#include <vector>

#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/quad4_rz.hpp"

namespace fuelsim {

struct M1Parameters final {
    double fuel_radius;
    double cladding_inner_radius;
    double cladding_outer_radius;
    double length;

    std::size_t fuel_radial_elements;
    std::size_t cladding_radial_elements;
    std::size_t axial_elements;

    ThermoelasticProperties fuel;
    ThermoelasticProperties cladding;
    double volumetric_heat_source;
    double outer_temperature;
    double initial_temperature;

    double gap_conductivity;
    double minimum_gap;
    double contact_penalty;
};

struct InterfaceSummary final {
    double minimum_gap;
    double maximum_gap;
    double maximum_contact_pressure;
    double total_heat_rate;
    double total_contact_force;
};

class M1Problem final : public NonlinearProblem {
  public:
    explicit M1Problem(M1Parameters parameters);

    const M1Parameters& parameters() const noexcept;
    const StructuredRzMesh& fuel_mesh() const noexcept;
    const StructuredRzMesh& cladding_mesh() const noexcept;
    const DofMap& dof_map() const noexcept;
    const Quad4RzThermoelasticKernel& fuel_kernel() const noexcept;
    const Quad4RzThermoelasticKernel& cladding_kernel() const noexcept;
    const Line2RzGapContactKernel& interface_kernel() const noexcept;

    std::size_t fuel_node_count() const noexcept;
    std::size_t cladding_node_offset() const noexcept;
    std::size_t fuel_global_node(std::size_t local_node) const;
    std::size_t cladding_global_node(std::size_t local_node) const;

    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    std::size_t fuel_element_count() const noexcept;
    std::size_t cladding_element_count() const noexcept;
    std::size_t interface_count() const noexcept;

    const std::vector<DirichletCondition>&
    dirichlet_conditions() const noexcept override;
    std::vector<double> initial_state() const;

    LocalDofs contribution_dofs(std::size_t contribution_index) const override;
    LocalResidual
    contribution_residual(std::size_t contribution_index,
                          const LocalValues& state) const override;
    LocalSystem linearize_contribution(std::size_t contribution_index,
                                       const LocalValues& state) const override;

    LocalDofs fuel_element_dofs(std::size_t element_index) const;
    LocalDofs cladding_element_dofs(std::size_t element_index) const;
    LocalDofs interface_dofs(std::size_t interface_index) const;

    const Line2RzInterfaceGeometry&
    interface_geometry(std::size_t interface_index) const;
    InterfaceSummary
    summarize_interface(const std::vector<double>& state) const;

  protected:
    void add_state_independent_residual(
        std::vector<double>& residual) const override;

  private:
    void build_geometries();
    void build_dirichlet_conditions();

    M1Parameters parameters_;
    StructuredRzMesh fuel_mesh_;
    StructuredRzMesh cladding_mesh_;
    DofMap dof_map_;
    Quad4RzThermoelasticKernel fuel_kernel_;
    Quad4RzThermoelasticKernel cladding_kernel_;
    Line2RzGapContactKernel interface_kernel_;
    std::vector<Quad4RzGeometry> fuel_geometries_;
    std::vector<Quad4RzGeometry> cladding_geometries_;
    std::vector<Line2RzInterfaceGeometry> interface_geometries_;
    std::vector<DirichletCondition> dirichlet_conditions_;
};

} // namespace fuelsim

#endif
