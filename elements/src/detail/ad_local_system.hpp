#pragma once
#include <adlite/adlite.hpp>
#include <cstddef>

namespace fuelsim::ad_local_system {
inline void make_passive(const double* values, std::size_t size, adlite::Scalar* result) {
    for (std::size_t entry = 0; entry < size; ++entry)
        result[entry] = values[entry];
}

inline void make_active(const double* values, std::size_t size, adlite::Scalar* result) {
    adlite::seed_identity(values, size, result);
}

inline void extract_residual(const adlite::Scalar* residual, std::size_t size, double* result) {
    for (std::size_t row = 0; row < size; ++row)
        result[row] = residual[row].value();
}

inline void
extract_system(const adlite::Scalar* residual, std::size_t size, double* residual_values, double* jacobian) {
    adlite::extract_jacobian(residual, size, size, residual_values, jacobian);
}
} // namespace fuelsim::ad_local_system
