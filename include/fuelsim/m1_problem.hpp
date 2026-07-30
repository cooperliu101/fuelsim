#ifndef FUELSIM_M1_PROBLEM_HPP
#define FUELSIM_M1_PROBLEM_HPP

#include <array>
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
    double fuel_length;
    double cladding_length;

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

struct ContactNodeSummary final {
    double z;
    bool projected;
    double gap;
    double pressure;
    double tributary_area;
    double tributary_length;
    double contact_force;
};

struct InterfaceSummary final {
    double minimum_gap;
    double maximum_gap;
    double minimum_contact_gap;
    double maximum_contact_pressure;
    double total_heat_rate;
    double total_contact_force;
    std::size_t projected_contact_nodes;
    std::size_t active_contact_nodes;
    double active_contact_length;
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
    const Line2RzGapHeatKernel& gap_heat_kernel() const noexcept;
    const NodeToLineRzContactKernel& contact_kernel() const noexcept;

    void set_volumetric_heat_source(double volumetric_heat_source);

    std::size_t fuel_node_count() const noexcept;
    std::size_t cladding_node_offset() const noexcept;
    std::size_t fuel_global_node(std::size_t local_node) const;
    std::size_t cladding_global_node(std::size_t local_node) const;

    std::size_t dof_count() const noexcept override;
    std::size_t contribution_count() const noexcept override;
    std::size_t fuel_element_count() const noexcept;
    std::size_t cladding_element_count() const noexcept;
    std::size_t thermal_interface_count() const noexcept;
    std::size_t contact_contribution_count() const noexcept;

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
    LocalDofs thermal_interface_dofs(std::size_t interface_index) const;
    LocalDofs contact_contribution_dofs(std::size_t contact_index) const;

    const Line2RzHeatGeometry&
    thermal_interface_geometry(std::size_t interface_index) const;
    const NodeToLineRzContactGeometry&
    contact_contribution_geometry(std::size_t contact_index) const;

    std::vector<ContactNodeSummary>
    summarize_contact_nodes(const std::vector<double>& state) const;
    InterfaceSummary
    summarize_interface(const std::vector<double>& state) const;

  protected:
    void add_state_independent_residual(
        std::vector<double>& residual) const override;

  private:
    void build_geometries();
    void build_dirichlet_conditions();

    M1Parameters _parameters;
    StructuredRzMesh _fuel_mesh;
    StructuredRzMesh _cladding_mesh;
    DofMap _dof_map;
    Quad4RzThermoelasticKernel _fuel_kernel;
    Quad4RzThermoelasticKernel _cladding_kernel;
    Line2RzGapHeatKernel _gap_heat_kernel;
    NodeToLineRzContactKernel _contact_kernel;
    std::vector<Quad4RzGeometry> _fuel_geometries;
    std::vector<Quad4RzGeometry> _cladding_geometries;
    std::vector<std::array<std::size_t, 4>> _thermal_interface_nodes;
    std::vector<Line2RzHeatGeometry> _thermal_interface_geometries;
    std::vector<std::array<std::size_t, 4>> _contact_contribution_nodes;
    std::vector<NodeToLineRzContactGeometry> _contact_geometries;
    std::vector<std::size_t> _contact_secondary_axial_indices;
    std::vector<DirichletCondition> _dirichlet_conditions;
};

} // namespace fuelsim

#endif
