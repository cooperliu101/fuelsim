#include "solver/petsc_solver.hpp"
#include "solver/section_modes.hpp"
#include "support/section_fixture.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
using namespace fuelsim;
using fuelsim::test::rectangle;
using fuelsim::test::region;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void verify_homogeneous() {
    constexpr double young = 70.0e9, h = 0.002, width = 0.02;
    const auto section = rectangle({-h / 2.0, h / 2.0}, width, {region("al", young, 0.3)});
    const auto basis = build_classic_section_basis(section);
    const auto response = evaluate_classic_section(section, basis, {0.01, 0.02, -0.03});
    const std::array<double, 3> expected{young * h * width,
        young * width * h * h * h / 12.0,
        young * h * width * width * width / 12.0};
    for (std::size_t i = 0; i < 3; ++i) {
        const double error = std::abs(response.tangent[4 * i] / expected[i] - 1.0);
        std::cout << "classic_stiffness_" << i << "=" << response.tangent[4 * i] << " relative_error=" << error
                  << " auxiliary_residual=" << basis.modes[i].auxiliary_relative_residual << '\n';
        require(error < 1.0e-10, "Uniform EA/EI failed");
        for (std::size_t j = 0; j < 3; ++j)
            if (i != j)
                require(std::abs(response.tangent[3 * i + j]) / std::sqrt(expected[i] * expected[j]) < 1.0e-10,
                    "Centered symmetric section must decouple extension and bending");
    }
    const auto rotation = classic_section_kinematics(section, basis, 0, {0.1, 0.2, -0.3}, {0.0, 0.4, -0.2}, {}, {});
    require(rotation.gradient[2] == 0.4 && rotation.gradient[6] == -0.4,
        "Thickness bending rotation must cancel axial warping gradient");
    require(rotation.gradient[5] == -0.2 && rotation.gradient[7] == 0.2,
        "Width bending rotation must cancel axial warping gradient");
    double stress_error = 0.0;
    for (std::size_t q = 0; q < section.points().size(); ++q) {
        const auto& p = section.points()[q];
        const double axial = 0.01 - p.position.x * 0.02 + p.position.y * 0.03;
        stress_error = std::max(stress_error, std::abs(response.stress[q][2] / (young * axial) - 1.0));
        require(std::abs(response.strain[q][2] - axial) < 1.0e-15, "Axial curvature sign failed");
        require(response.strain[q][4] == 0.0 && response.strain[q][5] == 0.0, "Classic shear cancellation failed");
        require(std::abs(response.stress[q][0]) < 0.1 && std::abs(response.stress[q][1]) < 0.1,
            "Poisson relaxation must eliminate transverse stress in this homogeneous patch");
    }
    std::cout << "classic_axial_stress_relative_error=" << stress_error << '\n';
    require(stress_error < 1.0e-10, "Axial stress patch failed");
    // A changing extension strain produces transverse correction shear; it must not be dropped.
    const auto varying = classic_section_strain(section, basis, 0, {}, {1.0, 0.0, 0.0});
    require(std::abs(varying[5] + 0.15 * section.points()[0].position.x) < 1.0e-12,
        "Axial derivative of Poisson contraction is missing");
    const std::array<double, 3> direction{0.3, -0.4, 0.2};
    std::array<double, 3> plus{0.01, 0.02, -0.03}, minus = plus;
    constexpr double step = 1.0e-6;
    for (std::size_t i = 0; i < 3; ++i) {
        plus[i] += step * direction[i];
        minus[i] -= step * direction[i];
    }
    const auto rp = evaluate_classic_section(section, basis, plus);
    const auto rm = evaluate_classic_section(section, basis, minus);
    for (std::size_t i = 0; i < 3; ++i) {
        double exact = 0.0;
        for (std::size_t j = 0; j < 3; ++j)
            exact += response.tangent[3 * i + j] * direction[j];
        require(std::abs((rp.force[i] - rm.force[i]) / (2.0 * step) - exact) < 1.0e-6 * expected[i],
            "Section tangent directional derivative failed");
    }
    double work = 0.0;
    const std::array<double, 3> amplitude{0.01, 0.02, -0.03};
    for (std::size_t i = 0; i < 3; ++i)
        work += 0.5 * amplitude[i] * response.force[i];
    require(std::abs(work / response.energy - 1.0) < 1.0e-12, "Section energy consistency failed");
}

void verify_gauge(const CrossSection& section, const ClassicSectionBasis& basis) {
    const auto n = section.nodes().size();
    for (const auto& mode : basis.modes) {
        double mx = 0.0, my = 0.0, mr = 0.0, scale = 0.0;
        for (const auto& point : section.points()) {
            double ux = 0.0, uy = 0.0;
            for (std::size_t i = 0; i < 8; ++i) {
                ux += point.shape[i] * mode.correction[point.nodes[i]];
                uy += point.shape[i] * mode.correction[n + point.nodes[i]];
            }
            const double x = point.position.x - section.centroid().x;
            const double y = point.position.y - section.centroid().y;
            mx += point.weight * ux;
            my += point.weight * uy;
            mr += point.weight * (-y * ux + x * uy);
            scale += point.weight * (std::abs(ux) + std::abs(uy));
        }
        require(std::abs(mx) + std::abs(my) + std::abs(mr) < 1.0e-10 * scale + 1.0e-25,
            "Poisson correction must be orthogonal to all three rigid section motions");
    }
}

void verify_invalid_inputs() {
    bool rejected = false;
    try {
        CrossSection empty({}, {}, {});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "Empty section must be rejected");
    rejected = false;
    try {
        (void)rectangle({0.001, -0.001}, 0.02, {region("inverted", 70e9, 0.3)});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "Inverted section must be rejected");
    auto section = rectangle({-0.001, 0.001}, 0.02, {region("valid", 70e9, 0.3)});
    auto basis = build_classic_section_basis(section);
    verify_gauge(section, basis);
    basis.modes[0].correction.clear();
    rejected = false;
    try {
        (void)evaluate_classic_section(section, basis, {});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "Mismatched basis must be rejected");
}

void verify_sandwich() {
    const std::vector<double> x{-0.002, -0.001, 0.001, 0.003};
    const std::array<double, 3> young{70e9, 90e9, 70e9};
    constexpr double width = 0.02;
    // Zero Poisson ratio gives an exact, independent heterogeneous reference without interface incompatibility.
    const auto section = rectangle(x,
        width,
        {region("al_bottom", young[0], 0.0), region("fuel", young[1], 0.0), region("al_top", young[2], 0.0)});
    double ea = 0.0, first = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        ea += young[i] * width * (x[i + 1] - x[i]);
        first += young[i] * width * (x[i + 1] * x[i + 1] - x[i] * x[i]) / 2.0;
    }
    const double center = first / ea;
    double ei = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
        ei += young[i] * width * (std::pow(x[i + 1] - center, 3) - std::pow(x[i] - center, 3)) / 3.0;
    const auto basis = build_classic_section_basis(section);
    const auto response = evaluate_classic_section(section, basis, {0.001, 1.0, 0.0});
    require(std::abs(section.elastic_center().x - center) < 1.0e-15, "Sandwich elastic center failed");
    require(std::abs(response.tangent[0] / ea - 1.0) < 1.0e-12, "Sandwich EA failed");
    require(std::abs(response.tangent[4] / ei - 1.0) < 1.0e-12, "Sandwich EI failed");
    for (std::size_t q = 0; q < section.points().size(); ++q) {
        const auto& p = section.points()[q];
        const double expected = young[p.region] * (0.001 - p.position.x + center);
        require(std::abs(response.stress[q][2] - expected) < 1.0e-6, "Sandwich regional stress failed");
    }
    std::cout << "sandwich_EA=" << ea << " EI=" << ei << " neutral_axis=" << center << '\n';
    const auto mixed = rectangle(x,
        width,
        {region("al_bottom", young[0], 0.33), region("fuel", young[1], 0.25), region("al_top", young[2], 0.33)},
        false,
        3);
    const auto mixed_basis = build_classic_section_basis(mixed);
    verify_gauge(mixed, mixed_basis);
    const auto mixed_response = evaluate_classic_section(mixed, mixed_basis, {0.001, 1.0, 0.0});
    require(mixed_response.energy > 0.0, "Mixed Poisson section energy must be positive");
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            require(std::abs(mixed_response.tangent[3 * i + j] - mixed_response.tangent[3 * j + i])
                        < 1.0e-10 * std::sqrt(mixed_response.tangent[4 * i] * mixed_response.tangent[4 * j]),
                "Mixed material section tangent must be symmetric");
    std::cout << "mixed_poisson_EA=" << mixed_response.tangent[0] << " EI=" << mixed_response.tangent[4]
              << " coupling=" << mixed_response.tangent[1] << '\n';
    const auto skew = rectangle({-0.001, 0.001}, width, {region("skew", 70e9, 0.3)}, true);
    const auto skew_basis = build_classic_section_basis(skew);
    const auto skew_response = evaluate_classic_section(skew, skew_basis, {1.0, 0.0, 0.0});
    require(std::abs(skew_response.tangent[0] / (70e9 * 0.002 * width) - 1.0) < 1.0e-10,
        "Nonrectangular affine section patch failed");
}
} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "Classic section internal contracts");
        std::cout << std::setprecision(16);
        verify_homogeneous();
        verify_sandwich();
        verify_invalid_inputs();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
