#pragma once
#include "material_functions.hpp"
#include <array>
#include <vector>

namespace fuelsim {
struct AxisymmetricStress final {
    adlite::Scalar rr, zz, hoop, rz;
};

struct AxisymmetricStressValues final {
    double rr, zz, hoop, rz;
};

struct AxisymmetricRotation final {
    adlite::Scalar rr{1.0}, rz{0.0}, zr{0.0}, zz{1.0}, hoop{1.0};
};

struct MaterialPointState final {
    std::array<double, 4> elastic_strain{}, plastic_strain{}, creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0;
    AxisymmetricStressValues stress{};
};

struct MaterialPointTrialState final {
    std::array<adlite::Scalar, 4> elastic_strain{}, plastic_strain{}, creep_strain{};
    adlite::Scalar equivalent_plastic_strain{0.0}, equivalent_creep_strain{0.0};
};

struct InelasticStressResponse final {
    AxisymmetricStress stress;
    MaterialPointTrialState trial_state;
};
} // namespace fuelsim

namespace fuelsim {
struct SymmetricTensor3Values final {
    double xx, yy, zz, xy, yz, xz;
};

struct CartesianMaterialPointState final {
    std::array<double, 6> elastic_strain{}, plastic_strain{}, creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0;
    SymmetricTensor3Values stress{};
};

struct CartesianInelasticStressResponse final {
    SymmetricTensor3 stress;
    CartesianMaterialPointState trial_state;
};

struct CartesianRotation final {
    adlite::Scalar xx{1.0}, xy{0.0}, xz{0.0};
    adlite::Scalar yx{0.0}, yy{1.0}, yz{0.0};
    adlite::Scalar zx{0.0}, zy{0.0}, zz{1.0};
};

using CartesianMaterialHistory = std::vector<CartesianMaterialPointState>;
} // namespace fuelsim
