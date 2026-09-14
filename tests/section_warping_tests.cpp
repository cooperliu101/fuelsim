#include "solver/petsc_solver.hpp"
#include "solver/section_warping.hpp"
#include "support/section_fixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
using namespace fuelsim;
using fuelsim::test::rectangle;
using fuelsim::test::region;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

double rectangle_torsion(double thickness, double width) {
    // NPTEL End Torsion of Prismatic Bars, equation 9.92, converted to full side lengths.
    const double pi = std::acos(-1.0);
    double sum = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const double odd = static_cast<double>(2 * i + 1);
        sum += std::tanh(odd * pi * width / (2.0 * thickness)) / std::pow(odd, 5);
    }
    return width * std::pow(thickness, 3) / 3.0 * (1.0 - 192.0 * thickness * sum / (std::pow(pi, 5) * width));
}

void verify_torsion() {
    constexpr double h = 0.002, width = 0.02, young = 70e9, nu = 0.3;
    const double reference = young / (2.0 * (1.0 + nu)) * rectangle_torsion(h, width);
    double previous = 1.0;
    for (int refinement : {1, 2, 4}) {
        std::vector<double> x;
        for (int i = 0; i <= refinement; ++i)
            x.push_back(-h / 2.0 + h * static_cast<double>(i) / static_cast<double>(refinement));
        const auto section = rectangle(x, width, {region("al", young, nu)}, false, 10 * refinement);
        const auto mode = build_section_torsion_mode(section);
        const double error = mode.warping.stiffness / reference - 1.0;
        std::cout << "torsion_refinement=" << refinement << " nodes=" << section.nodes().size()
                  << " GJ=" << mode.warping.stiffness << " reference_GJ=" << reference << " relative_error=" << error
                  << " residual=" << mode.warping.relative_residual << '\n';
        require(error > -1.0e-10 && error < previous, "Torsion primal energy must converge from above");
        previous = error;
        double mean = 0.0, torque = 0.0, work = 0.0;
        for (std::size_t i = 0; i < mode.transverse.size(); ++i)
            work += mode.transverse[i] * mode.warping.condensed_force[i];
        for (std::size_t q = 0; q < section.points().size(); ++q) {
            const auto& point = section.points()[q];
            const auto k = section_torsion_kinematics(section, mode, q, 0.0, 1.0, 2.0);
            const double gx = k.gradient[2] + k.gradient[6], gy = k.gradient[5] + k.gradient[7];
            mean += point.weight * k.displacement[2];
            torque += point.weight * young / (2.0 * (1.0 + nu)) * (point.position.x * gy - point.position.y * gx);
            require(std::abs(k.gradient[8] - 2.0 * k.displacement[2]) < 1.0e-15,
                "Nonuniform torsion must retain axial warping strain");
        }
        require(std::abs(mean) < section.area() * 1.0e-15, "Warping must have zero integral mean");
        require(std::abs(torque / mode.warping.stiffness - 1.0) < 1.0e-9, "Stress torque must equal energy stiffness");
        require(std::abs(work / mode.warping.stiffness - 1.0) < 1.0e-9, "Condensed operator energy failed");
    }
    require(previous < 1.0e-4, "Refined rectangle GJ must agree within 0.01 percent");
}

void verify_condensation() {
    const auto section = rectangle({-0.001, 0.001}, 0.02, {region("al", 70e9, 0.3)});
    SectionWarpingSolver solver(section);
    const auto n = section.nodes().size();
    std::vector<double> translation(2 * n), dilation(2 * n), rotation(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        translation[i] = 1.0;
        dilation[i] = section.nodes()[i].x;
        dilation[n + i] = section.nodes()[i].y;
        rotation[i] = -section.nodes()[i].y;
        rotation[n + i] = section.nodes()[i].x;
    }
    const auto shift = solver.solve(translation);
    const auto gradient = solver.solve(dilation);
    const auto twist = solver.solve(rotation);
    require(shift.stiffness < 1.0e-15, "Translation must be cancellable by linear warping");
    require(gradient.stiffness < 1.0e-20, "Gradient transverse field must expose nonrigid K_S kernel");
    const double step = 1.0e-5;
    std::vector<double> plus = rotation, minus = rotation;
    for (std::size_t i = 0; i < plus.size(); ++i) {
        plus[i] += step * translation[i];
        minus[i] -= step * translation[i];
    }
    const auto rp = solver.solve(plus), rm = solver.solve(minus);
    double expected = 0.0;
    for (std::size_t i = 0; i < plus.size(); ++i)
        expected += twist.condensed_force[i] * translation[i];
    require(std::abs((rp.stiffness - rm.stiffness) / (4.0 * step) - expected) < 1.0e-7,
        "Condensed force must differentiate minimized shear energy");
    bool rejected = false;
    try {
        (void)solver.solve({});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "Invalid transverse field size must be rejected");
}
} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "Section warping internal contracts");
        std::cout << std::setprecision(16);
        verify_torsion();
        verify_condensation();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
