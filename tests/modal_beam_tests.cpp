#include "solver/petsc_solver.hpp"
#include "solver/reduced_section_basis.hpp"
#include "support/section_fixture.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void verify() {
    const auto section = fuelsim::test::rectangle({-0.001, 0.001}, 0.02, {fuelsim::test::region("al", 70e9, 0.3)});
    const auto basis = fuelsim::build_reduced_section_basis(section, true, 2);
    const auto element = fuelsim::make_modal_beam_element(0.0, 0.5);
    const auto n = 6 * basis.modes.size();
    std::vector<double> state(n);
    constexpr double epsilon = 0.001;
    state[1] = epsilon;
    state[3] = 0.5 * epsilon;
    state[4] = epsilon;
    const auto extension = fuelsim::evaluate_modal_beam(section, basis, element, state);
    const double extension_energy = 0.5 * 2.8e6 * epsilon * epsilon * 0.5;
    std::cout << "beam_extension_energy_error=" << std::abs(extension.energy / extension_energy - 1.0) << '\n';
    require(std::abs(extension.energy / extension_energy - 1.0) < 1.0e-10, "Beam extension energy failed");
    state.assign(n, 0.0);
    state[8] = 0.02;
    state[9] = 0.02 * 0.5 * 0.5 / 2.0;
    state[10] = 0.02 * 0.5;
    state[11] = 0.02;
    const auto bending = fuelsim::evaluate_modal_beam(section, basis, element, state);
    const double bending_energy = 0.5 * (0.9333333333333333) * 0.02 * 0.02 * 0.5;
    std::cout << "beam_bending_energy_error=" << std::abs(bending.energy / bending_energy - 1.0) << '\n';
    require(std::abs(bending.energy / bending_energy - 1.0) < 1.0e-10, "Beam bending energy failed");
    // All six three-dimensional rigid motions must remain zero-energy modes.
    for (std::size_t rigid = 0; rigid < 6; ++rigid) {
        state.assign(n, 0.0);
        if (rigid < 3) {
            state[6 * rigid] = state[6 * rigid + 3] = 1.0;
        } else if (rigid < 5) {
            const auto offset = 6 * (rigid - 2);
            state[offset + 1] = state[offset + 4] = 1.0;
            state[offset + 3] = 0.5;
        } else {
            state[18] = state[21] = 1.0;
        }
        const auto response = fuelsim::evaluate_modal_beam(section, basis, element, state, false);
        require(response.energy < 1.0e-16, "Beam rigid motion acquired strain energy");
    }
    for (std::size_t i = 0; i < n; ++i)
        state[i] = 1.0e-4 * std::sin(static_cast<double>(i + 1));
    const auto response = fuelsim::evaluate_modal_beam(section, basis, element, state);
    const auto residual_only = fuelsim::evaluate_modal_beam(section, basis, element, state, false);
    require(response.residual == residual_only.residual, "Beam residual-only and tangent paths disagree");
    std::vector<double> direction(n), plus = state, minus = state;
    constexpr double step = 1.0e-7;
    for (std::size_t i = 0; i < n; ++i) {
        direction[i] = std::cos(static_cast<double>(i + 1));
        plus[i] += step * direction[i];
        minus[i] -= step * direction[i];
    }
    const auto rp = fuelsim::evaluate_modal_beam(section, basis, element, plus, false);
    const auto rm = fuelsim::evaluate_modal_beam(section, basis, element, minus, false);
    double error = 0.0, norm = 0.0, work = 0.0, symmetry = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double exact = 0.0;
        work += 0.5 * state[i] * response.residual[i];
        for (std::size_t j = 0; j < n; ++j) {
            exact += response.jacobian[n * i + j] * direction[j];
            symmetry = std::max(symmetry,
                std::abs(response.jacobian[n * i + j] - response.jacobian[n * j + i])
                    / std::sqrt(response.jacobian[n * i + i] * response.jacobian[n * j + j]));
        }
        error = std::hypot(error, (rp.residual[i] - rm.residual[i]) / (2.0 * step) - exact);
        norm = std::hypot(norm, exact);
    }
    std::cout << "beam_directional_derivative_error=" << error / norm << " symmetry=" << symmetry << '\n';
    require(error / norm < 1.0e-8 && symmetry < 1.0e-10, "Beam consistent tangent failed");
    require(std::abs(work / response.energy - 1.0) < 1.0e-10, "Beam energy and residual virtual work disagree");
    require(response.strain.size() == 6 * section.points().size(), "Section material points were reduced away");
}
} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "Modal beam local contracts");
        verify();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
