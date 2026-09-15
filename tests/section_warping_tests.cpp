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
    // Curl of (a^2-x^2)*(b^2-y^2): divergence-free and tangent to all four
    // boundaries, so its weak warping load vanishes although its shear energy
    // is nonzero. This checks Neumann compatibility at a nearly zero RHS.
    std::vector<double> solenoidal(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = section.nodes()[i];
        solenoidal[i] = -2.0 * p.y * (1e-6 - p.x * p.x);
        solenoidal[n + i] = 2.0 * p.x * (1e-4 - p.y * p.y);
    }
    const auto closed = solver.solve(solenoidal);
    for (double value : closed.axial)
        require(std::abs(value) < 1e-20, "Solenoidal transverse field acquired nonzero warping");
    require(closed.relative_residual < 1e-9 && closed.stiffness > 0.0,
        "Near-zero Neumann load lost equilibrium or physical shear energy");
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

void verify_axial_corrector() {
    std::vector<double> x;
    for (int i = 0; i <= 8; ++i)
        x.push_back(-0.001 + 0.002 * static_cast<double>(i) / 8.0);
    const auto section = rectangle(x, 0.02, {region("al", 70e9, 0.3)});
    SectionWarpingSolver solver(section);
    std::vector<double> source;
    for (const auto& node : section.nodes())
        source.push_back(node.x);
    const auto result = solver.solve_axial_corrector(source);
    double difference = 0.0, norm = 0.0, mean = 0.0;
    for (const auto& point : section.points()) {
        double actual = 0.0;
        for (std::size_t i = 0; i < 8; ++i)
            actual += point.shape[i] * result.axial[point.nodes[i]];
        const double coordinate = point.position.x;
        // -G*chi,xx=Czz*x and chi,x=0 at x=+-a. Czz/G=3.5 for nu=0.3.
        const double expected = 3.5 * (1e-6 * coordinate / 2.0 - coordinate * coordinate * coordinate / 6.0);
        difference += point.weight * (actual - expected) * (actual - expected);
        norm += point.weight * expected * expected;
        mean += point.weight * actual;
    }
    std::cout << "axial_corrector_L2=" << std::sqrt(difference / norm) << " residual=" << result.relative_residual
              << '\n';
    require(std::sqrt(difference / norm) < 0.002, "Axial corrector failed the independent Neumann polynomial solution");
    require(std::abs(mean) < 1e-24, "Axial corrector gauge is not zero mean");
    const auto constant = solver.solve_axial_corrector(std::vector<double>(source.size(), 123.0));
    for (double value : constant.axial)
        require(value == 0.0, "A constant axial source must be removed by the Neumann compatibility projection");
}
} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "Section warping internal contracts");
        std::cout << std::setprecision(16);
        verify_torsion();
        verify_condensation();
        verify_axial_corrector();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
