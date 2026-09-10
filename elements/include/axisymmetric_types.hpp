#pragma once
#include "material_functions.hpp"
#include <array>

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
