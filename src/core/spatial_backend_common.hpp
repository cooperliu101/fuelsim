#pragma once
#include "spatial_backends.hpp"

namespace fuelsim {
std::string set_commit_failure(std::vector<double>& values, const std::exception_ptr& failure);
void sum_commit_partitions(std::vector<double>& values,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    std::string failure_message);
void synchronize_commit_diagnostics(std::vector<double>& residual,
    std::vector<double>& external,
    TransientConservationSummary& conservation,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    const std::exception_ptr& failure);
void synchronize_cartesian_commit(std::vector<std::vector<CartesianMaterialHistory>>& histories,
    std::vector<double>& residual,
    std::vector<double>& external,
    TransientConservationSummary& conservation,
    std::size_t first,
    std::size_t last,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    const std::exception_ptr& failure);
void synchronize_rz_commit(const std::vector<MaterialPointState*>& points,
    std::size_t points_per_element,
    std::vector<double>& residual,
    std::vector<double>& external,
    TransientConservationSummary& conservation,
    std::size_t first,
    std::size_t last,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    const std::exception_ptr& failure);
double body_point_work(const RegionDefinition& region,
    const MaterialFunctionContext& context,
    double reference_measure,
    const std::array<double, 3>& displacement_increment);
void finalize_conservation(const NonlinearProblem& problem,
    const std::vector<double>& current,
    const std::vector<double>& old,
    const std::vector<double>& current_raw_residual,
    const std::vector<double>& old_raw_residual,
    const std::vector<double>& current_external_load_residual,
    const std::vector<double>& old_external_load_residual,
    TransientConservationSummary& result);
Cax4LocalValues rz_local_values(const std::vector<double>& values);
Cax4LocalValues
gather_rz_state(const rz::SpatialAssembly& spatial, std::size_t index, const std::vector<double>& global_state);
bool finite_stress(const AxisymmetricStressValues& stress);
bool valid_material_state(const MaterialPointState& state);
bool valid_material_state(const CartesianMaterialPointState& state);

namespace rz {
double stress_strain_inner_product(const AxisymmetricStressValues& stress,
    const std::array<double, 4>& strain) noexcept;
void accumulate_material_conservation(TransientConservationSummary& result,
    const MaterialPointState& old,
    const MaterialPointState& current,
    double measure);
} // namespace rz

namespace cartesian {
double stress_strain_inner_product(const SymmetricTensor3Values& stress, const std::array<double, 6>& strain) noexcept;
std::array<double, 6> strain_difference(const std::array<double, 6>& current, const std::array<double, 6>& old);
double trapezoidal_stress_strain_inner_product(const SymmetricTensor3Values& old_stress,
    const SymmetricTensor3Values& new_stress,
    const std::array<double, 6>& strain_increment) noexcept;
SymmetricTensor3Values rotate_tensor_values(const std::array<double, 6>& tensor, const CartesianRotation& rotation);
SymmetricTensor3Values rotate_tensor_values(const SymmetricTensor3Values& tensor, const CartesianRotation& rotation);
std::array<double, 6> components(const SymmetricTensor3Values& tensor);
} // namespace cartesian
} // namespace fuelsim
