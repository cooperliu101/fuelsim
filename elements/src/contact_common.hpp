#pragma once
#include "contact_types.hpp"

namespace fuelsim::contact_common {
using ActivePoint3 = std::array<adlite::Scalar, 3>;

struct FrictionResult final {
    bool sliding = false;
    ActivePoint3 tangential_slip{}, elastic_tangential_slip{}, tangential_traction_vector{};
    adlite::Scalar tangential_traction{0.0}, tangential_force{0.0}, friction_dissipation{0.0};
};

FrictionResult friction_return(const NormalContactProperties& properties,
    const ContactPointHistory& history,
    const ActivePoint3& transported_history,
    const ActivePoint3& transported_total_history,
    const ActivePoint3& relative_increment,
    const ActivePoint3& normal,
    const adlite::Scalar& pressure,
    const adlite::Scalar& tributary_area);
// Nonfinite states throw; false means a finite Newton increment is still too large.
bool projection_increment_converged(const adlite::Scalar& delta_xi,
    const adlite::Scalar& delta_eta,
    const adlite::Scalar& xi,
    const adlite::Scalar& eta,
    std::size_t maximum_width);
// At either max-law corner ADlite supplies the mean of the two branch derivatives.
adlite::Scalar gap_conductance(const GapHeatProperties& properties,
    const adlite::Scalar& gap,
    const adlite::Scalar& secondary_temperature,
    const adlite::Scalar& primary_temperature);
} // namespace fuelsim::contact_common
