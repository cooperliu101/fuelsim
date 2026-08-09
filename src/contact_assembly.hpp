#ifndef FUELSIM_CONTACT_ASSEMBLY_HPP
#define FUELSIM_CONTACT_ASSEMBLY_HPP

#include "fuelsim/interface.hpp"
#include "spatial_layout.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {

class SpatialAssembly;

class ContactAssembly final {
  public:
    ContactAssembly() = default;

    std::size_t contact_count(const SpatialLayout& layout) const noexcept;
    const ContactDefinition& contact(std::size_t contact_index,
                                     const SpatialLayout& layout) const;
    const std::vector<std::vector<ContactPointHistory>>&
    committed_histories() const noexcept;
    bool uses_augmented_contact(const SpatialLayout& layout) const noexcept;
    AugmentedContactUpdate
    update_augmented_multipliers(SpatialAssembly& assembly,
                                 const std::vector<double>& state,
                                 std::size_t completed_updates);
    void commit_state(SpatialAssembly& assembly,
                      const std::vector<double>& state);
    void restore_state(const SpatialAssembly& assembly,
                       const std::vector<double>& state,
                       std::vector<std::vector<ContactPointHistory>> histories);

  private:
    friend class SpatialAssembly;

    struct ThermalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        Line2RzHeatGeometry geometry;
    };
    struct MechanicalContribution final {
        std::size_t contact;
        std::array<std::size_t, 4> nodes;
        NodeToLineRzContactGeometry geometry;
        std::size_t secondary;
        std::size_t primary;
        mutable bool active = false;
    };

    std::vector<Line2RzGapHeatKernel> _thermal_kernels;
    std::vector<NodeToLineRzContactKernel> _mechanical_kernels;
    std::vector<ThermalContribution> _thermal_contributions;
    std::vector<MechanicalContribution> _mechanical_contributions;
    mutable std::vector<std::vector<bool>> _projected_mechanical_nodes;
    std::vector<std::vector<ContactPointHistory>> _contact_histories;
    std::vector<double> _committed_contact_solution;
    std::vector<SpatialLayout::ResolvedBoundary> _primary_boundaries;
    std::vector<SpatialLayout::ResolvedBoundary> _secondary_boundaries;
};

} // namespace fuelsim

#endif
