#ifndef FUELSIM_LOCAL_SYSTEM_HPP
#define FUELSIM_LOCAL_SYSTEM_HPP

#include <adlite/adlite.hpp>

#include <array>
#include <cstddef>

namespace fuelsim {

constexpr std::size_t local_dof_count = 12;
constexpr std::size_t local_jacobian_size = local_dof_count * local_dof_count;

using LocalDofs = std::array<std::size_t, local_dof_count>;
using LocalValues = std::array<double, local_dof_count>;
using LocalResidual = std::array<double, local_dof_count>;
using LocalJacobian = std::array<double, local_jacobian_size>;
using LocalAdValues = std::array<adlite::Scalar, local_dof_count>;

struct LocalSystem final {
    LocalResidual residual;
    LocalJacobian jacobian;
};

} // namespace fuelsim

#endif
