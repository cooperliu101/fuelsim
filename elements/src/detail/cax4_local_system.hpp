#pragma once
#include "ad_local_system.hpp"
#include "cax4_types.hpp"

namespace fuelsim::cax4_detail {
inline Cax4LocalAdValues ad_state(const Cax4LocalValues& state, bool active = false) {
    Cax4LocalAdValues result{};
    if (active)
        ad_local_system::make_active(state.data(), state.size(), result.data());
    else
        ad_local_system::make_passive(state.data(), state.size(), result.data());
    return result;
}

inline Cax4LocalResidual
values(const Cax4LocalAdValues& state, const Cax4LocalAdValues& residual, Cax4LocalJacobian* jacobian = nullptr) {
    Cax4LocalResidual result{};
    if (jacobian == nullptr)
        ad_local_system::extract_residual(residual.data(), residual.size(), result.data());
    else
        ad_local_system::extract_system(residual.data(), state.size(), result.data(), jacobian->data());
    return result;
}

} // namespace fuelsim::cax4_detail
