#pragma once
#include <adlite/adlite.hpp>
#include <array>

namespace fuelsim::cartesian_detail {
using Matrix3 = std::array<std::array<double, 3>, 3>;
using ActiveMatrix3 = std::array<std::array<adlite::Scalar, 3>, 3>;

double determinant(const Matrix3& matrix);
adlite::Scalar determinant(const ActiveMatrix3& matrix);
Matrix3 inverse(const Matrix3& matrix, double determinant_value);
ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& determinant_value);
ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second);

Matrix3 multiply(const Matrix3& first, const Matrix3& second);
} // namespace fuelsim::cartesian_detail
