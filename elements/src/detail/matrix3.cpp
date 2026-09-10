#include "matrix3.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::cartesian_detail {
double determinant(const Matrix3& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
           - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
           + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

adlite::Scalar determinant(const ActiveMatrix3& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
           - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
           + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

Matrix3 inverse(const Matrix3& matrix, double determinant_value) {
    return {{{{(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant_value,
                 (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant_value,
                 (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant_value}},
        {{(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant_value,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant_value,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant_value}},
        {{(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant_value,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant_value,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant_value}}}};
}

ActiveMatrix3 inverse(const ActiveMatrix3& matrix, const adlite::Scalar& determinant_value) {
    return {{{{(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant_value,
                 (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant_value,
                 (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant_value}},
        {{(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant_value,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant_value,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant_value}},
        {{(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant_value,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant_value,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant_value}}}};
}

ActiveMatrix3 multiply(const ActiveMatrix3& first, const Matrix3& second) {
    ActiveMatrix3 result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                result[i][j] += first[i][k] * second[k][j];
    return result;
}

Matrix3 multiply(const Matrix3& first, const Matrix3& second) {
    Matrix3 result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                result[i][j] += first[i][k] * second[k][j];
    return result;
}
} // namespace fuelsim::cartesian_detail
