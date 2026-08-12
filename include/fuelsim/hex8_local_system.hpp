#ifndef FUELSIM_HEX8_LOCAL_SYSTEM_HPP
#define FUELSIM_HEX8_LOCAL_SYSTEM_HPP

#include <adlite/adlite.hpp>

#include <array>
#include <cstddef>

namespace fuelsim {

inline constexpr std::size_t hex8_node_count = 8;
inline constexpr std::size_t hex8_local_dof_count = 32;
inline constexpr std::size_t quad4_face_node_count = 4;
inline constexpr std::size_t quad4_face_local_dof_count = 16;

using Hex8LocalDofs = std::array<std::size_t, hex8_local_dof_count>;
using Hex8LocalValues = std::array<double, hex8_local_dof_count>;
using Hex8LocalResidual = std::array<double, hex8_local_dof_count>;
using Hex8LocalJacobian = std::array<double, hex8_local_dof_count * hex8_local_dof_count>;
using Hex8LocalAdValues = std::array<adlite::Scalar, hex8_local_dof_count>;

struct Hex8LocalSystem final {
    Hex8LocalResidual residual;
    Hex8LocalJacobian jacobian;
};

using Quad4FaceLocalDofs = std::array<std::size_t, quad4_face_local_dof_count>;
using Quad4FaceLocalValues = std::array<double, quad4_face_local_dof_count>;
using Quad4FaceLocalResidual = std::array<double, quad4_face_local_dof_count>;
using Quad4FaceLocalJacobian = std::array<double, quad4_face_local_dof_count * quad4_face_local_dof_count>;
using Quad4FaceLocalAdValues = std::array<adlite::Scalar, quad4_face_local_dof_count>;

struct Quad4FaceLocalSystem final {
    Quad4FaceLocalResidual residual;
    Quad4FaceLocalJacobian jacobian;
};

} // namespace fuelsim

#endif
