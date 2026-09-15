#include "../src/solver/modal_condensation.hpp"
#include "../src/solver/modal_solid_end.hpp"
#include "solver/petsc_solver.hpp"
#include "solver/reduced_section_basis.hpp"
#include "solver/section_modes.hpp"
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

void verify_solid_interface() {
    const std::array<std::array<double, 3>, 20> locations{{{-1, -1, 0},
        {1, -1, 0},
        {1, 1, 0},
        {-1, 1, 0},
        {-1, -1, 1},
        {1, -1, 1},
        {1, 1, 1},
        {-1, 1, 1},
        {0, -1, 0},
        {1, 0, 0},
        {0, 1, 0},
        {-1, 0, 0},
        {-1, -1, 0.5},
        {1, -1, 0.5},
        {1, 1, 0.5},
        {-1, 1, 0.5},
        {0, -1, 1},
        {1, 0, 1},
        {0, 1, 1},
        {-1, 0, 1}}};
    fuelsim::Hex20Coordinates coordinates{};
    std::array<fuelsim::SolidDisplacementRow, 60> rows;
    std::size_t count = 18;
    for (std::size_t i = 0; i < 20; ++i) {
        const double x = 0.001 * locations[i][0], y = 0.01 * locations[i][1], z = 0.005 * locations[i][2];
        coordinates[i] = {x, y, z};
        const std::array<double, 6> trace{1.0, x, y, x * x, x * y, y * y};
        for (std::size_t c = 0; c < 3; ++c) {
            if (z == 0.0) {
                for (std::size_t j = 0; j < trace.size(); ++j)
                    if (trace[j] != 0.0)
                        rows[20 * c + i].emplace_back(6 * c + j, trace[j]);
            } else
                rows[20 * c + i].emplace_back(count++, 1.0);
        }
    }
    const auto element = fuelsim::make_modal_solid_element(coordinates, rows, 0, 0);
    const auto material = fuelsim::test::region("al", 70e9, 0.3);
    require(element.dofs.size() == count, "Solid trace lost an independent interface or interior unknown");
    for (int kind = 0; kind < 8; ++kind) {
        std::vector<double> state(count);
        // Six rigid motions, uniform extension, and exact shear-free bending.
        if (kind < 3)
            state[6 * static_cast<std::size_t>(kind)] = 1.0;
        if (kind == 3)
            state[14] = 1.0;
        if (kind == 4)
            state[13] = -1.0;
        if (kind == 5) {
            state[2] = -1.0;
            state[7] = 1.0;
        }
        if (kind == 6) {
            state[1] = -0.3;
            state[8] = -0.3;
        }
        if (kind == 7) {
            state[3] = 0.15;
            state[5] = -0.15;
            state[10] = 0.3;
        }
        for (std::size_t i = 0; i < 20; ++i) {
            const auto& p = coordinates[i];
            std::array<double, 3> u{};
            if (kind < 3)
                u[static_cast<std::size_t>(kind)] = 1.0;
            if (kind == 3)
                u = {0.0, -p.z, p.y};
            if (kind == 4)
                u = {p.z, 0.0, -p.x};
            if (kind == 5)
                u = {-p.y, p.x, 0.0};
            if (kind == 6)
                u = {-0.3 * p.x, -0.3 * p.y, p.z};
            if (kind == 7)
                u = {0.5 * p.z * p.z + 0.15 * (p.x * p.x - p.y * p.y), 0.3 * p.x * p.y, -p.x * p.z};
            if (p.z != 0.0)
                for (std::size_t c = 0; c < 3; ++c)
                    state[rows[20 * c + i][0].first] = u[c];
        }
        const auto response = fuelsim::evaluate_modal_solid(element, material, state, false, true);
        if (kind < 6)
            require(response.energy < 1e-16, "Solid/modal interface penalized a rigid motion");
        else {
            const double expected = 0.5 * 70e9 * 0.002 * 0.02 * 0.005 * (kind == 6 ? 1.0 : 0.002 * 0.002 / 12.0);
            require(std::abs(response.energy / expected - 1) < 1e-12, "Solid interface failed the EA/EI patch energy");
            for (std::size_t q = 0; q < response.strain.size(); ++q) {
                const auto& p = element.geometry.mechanical_points[q].position;
                require(std::abs(response.strain[q][2] - (kind == 6 ? 1.0 : -p.x)) < 1e-12
                            && std::abs(response.strain[q][4]) < 1e-12 && std::abs(response.strain[q][5]) < 1e-12,
                    "Solid interface extension/bending strain or shear cancellation failed");
            }
        }
    }
    std::vector<double> state(count), direction(count), plus(count), minus(count);
    for (std::size_t i = 0; i < count; ++i) {
        state[i] = 1e-5 * std::sin(static_cast<double>(i + 1));
        direction[i] = std::cos(static_cast<double>(i + 1));
        plus[i] = state[i] + 1e-7 * direction[i];
        minus[i] = state[i] - 1e-7 * direction[i];
    }
    const auto response = fuelsim::evaluate_modal_solid(element, material, state, true, true);
    const auto residual = fuelsim::evaluate_modal_solid(element, material, state, false, false);
    require(residual.residual == response.residual, "Solid end residual changed when requesting a tangent");
    const auto rp = fuelsim::evaluate_modal_solid(element, material, plus, false, false);
    const auto rm = fuelsim::evaluate_modal_solid(element, material, minus, false, false);
    double error = 0.0, norm = 0.0, work = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        double exact = 0.0;
        work += 0.5 * state[i] * response.residual[i];
        for (std::size_t j = 0; j < count; ++j)
            exact += response.jacobian[count * i + j] * direction[j];
        error = std::hypot(error, (rp.residual[i] - rm.residual[i]) / 2e-7 - exact);
        norm = std::hypot(norm, exact);
    }
    std::cout << "solid_interface_directional_derivative_error=" << error / norm << '\n';
    require(error / norm < 1e-8 && std::abs(work / response.energy - 1) < 1e-12,
        "Solid end interface tangent or virtual work failed");
}

void verify_reflection_tangent() {
    for (bool skew : {false, true}) {
        const auto section =
            fuelsim::test::rectangle({-0.001, 0.001}, 0.02, {fuelsim::test::region("al", 70e9, 0.3)}, skew, 2);
        const auto basis = fuelsim::build_reduced_section_basis(section, false, 0);
        require(basis.reflected_nodes[0].empty() == skew && basis.reflected_nodes[1].empty() == skew,
            "Reflection optimization assumed a symmetry absent from the actual mesh");
        const auto tangent = fuelsim::linearize_modal_section(section, basis);
        auto unpartitioned = basis;
        unpartitioned.reflected_nodes = {};
        const auto dense = fuelsim::linearize_modal_section(section, unpartitioned);
        const auto element = fuelsim::make_modal_beam_element(0.0, 0.05);
        const auto actual_stiffness = fuelsim::modal_beam_linear_stiffness(tangent, element);
        const auto full_stiffness = fuelsim::modal_beam_linear_stiffness(dense, element);
        // Classical q' bending shear is exactly zero, so its section diagonal
        // is only roundoff. Use the physical nodal energy scales, whose diagonal
        // entries are positive, rather than divide by that numerical zero.
        constexpr std::size_t size = 18;
        for (std::size_t i = 0; i < size; ++i)
            for (std::size_t j = 0; j < size; ++j) {
                const double scale = std::sqrt(full_stiffness[size * i + i] * full_stiffness[size * j + j]);
                require(std::abs(actual_stiffness[size * i + j] - full_stiffness[size * i + j]) <= 1e-11 * scale,
                    "Reflection separation changed the fully integrated material tangent");
            }
        const double bending_coupling = tangent.tangent[12 * 6 + 10];
        if (skew) {
            const auto classical = fuelsim::evaluate_classic_section(section,
                fuelsim::build_classic_section_basis(section),
                {0.0, 1.0, 1.0});
            require(std::abs(bending_coupling / classical.tangent[5] - 1.0) < 1e-10,
                "Skew section lost its physical cross-bending stiffness");
        } else {
            require(bending_coupling == 0.0, "Opposite reflection sectors retained a spurious matrix coupling");
        }
    }
    const auto layers = fuelsim::test::rectangle({-0.001, 0.0, 0.001},
        0.02,
        {fuelsim::test::region("al", 70e9, 0.3), fuelsim::test::region("fuel", 140e9, 0.25)},
        false,
        2);
    const auto layered_basis = fuelsim::build_reduced_section_basis(layers, false, 0);
    require(layered_basis.reflected_nodes[0].empty() && !layered_basis.reflected_nodes[1].empty(),
        "Reflection optimization ignored the actual material regions");
}

void verify_width_enrichment() {
    const auto section = fuelsim::test::rectangle({-0.001, -0.0005, 0.0, 0.0005, 0.001},
        0.02,
        {fuelsim::test::region("al", 70e9, 0.3)},
        false,
        8);
    const auto original = fuelsim::build_reduced_section_basis(section, true, 8);
    const auto seed = fuelsim::make_independent_section_basis(section, fuelsim::make_mixed_bending_basis(original));
    const auto inner = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double result = 0.0;
        const auto n = section.nodes().size();
        for (const auto& point : section.points())
            for (std::size_t c = 0; c < 3; ++c) {
                double x = 0.0, y = 0.0;
                for (std::size_t i = 0; i < 8; ++i) {
                    x += point.shape[i] * a[c * n + point.nodes[i]];
                    y += point.shape[i] * b[c * n + point.nodes[i]];
                }
                result += point.weight * x * y;
            }
        return result;
    };
    const auto element = fuelsim::make_modal_beam_element(0.0, 0.05);
    fuelsim::ReducedSectionBasis previous;
    for (int strips : {1, 2, 4}) {
        std::vector<double> lines;
        for (int i = 0; i <= strips; ++i)
            lines.push_back(-0.01 + 0.02 * i / strips);
        const auto basis = fuelsim::enrich_section_width(section, seed, lines);
        std::cout << "width_strips=" << strips << " modes=" << basis.modes.size() << '\n';
        for (std::size_t i = 0; i < basis.modes.size(); ++i) {
            for (std::size_t order = 1; order < 3; ++order)
                require(std::all_of(basis.modes[i].coefficient[order].begin(),
                            basis.modes[i].coefficient[order].end(),
                            [](double value) { return value == 0.0; }),
                    "Width line amplitudes lost an unrepresented derivative field");
            const auto& a = basis.modes[i].coefficient[0];
            for (std::size_t j = 0; j < i; ++j) {
                const auto& b = basis.modes[j].coefficient[0];
                require(std::abs(inner(a, b)) < 1e-10 * std::sqrt(inner(a, a) * inner(b, b)),
                    "Width basis retained a dependent or nonorthogonal direction");
            }
        }
        for (const auto& mode : previous.modes) {
            auto remainder = mode.coefficient[0];
            for (const auto& direction : basis.modes) {
                const auto& field = direction.coefficient[0];
                const double coefficient = inner(mode.coefficient[0], field) / inner(field, field);
                for (std::size_t i = 0; i < field.size(); ++i)
                    remainder[i] -= coefficient * field[i];
            }
            require(inner(remainder, remainder) < 1e-20 * inner(mode.coefficient[0], mode.coefficient[0]),
                "Refining width lines lost the coarser conforming displacement space");
        }
        for (std::size_t classic = 0; classic < 4; ++classic) {
            std::vector<double> reference(6 * original.modes.size()), state(6 * basis.modes.size());
            for (std::size_t side = 0; side < 2; ++side) {
                const double z = side == 0 ? element.lower : element.upper;
                const std::array<double, 5> derivative =
                    classic == 0 || classic == 3 ? std::array<double, 5>{0.001 * z, 0.001, 0.0, 0.0, 0.0}
                                                 : std::array<double, 5>{0.01 * z * z, 0.02 * z, 0.02, 0.0, 0.0};
                for (std::size_t order = 0; order < 3; ++order) {
                    reference[6 * classic + 3 * side + order] = derivative[order];
                    std::vector<double> field(3 * section.nodes().size());
                    for (std::size_t k = 0; k < 3; ++k)
                        for (std::size_t i = 0; i < field.size(); ++i)
                            field[i] += original.modes[classic].coefficient[k][i] * derivative[order + k];
                    for (std::size_t m = 0; m < basis.modes.size(); ++m) {
                        const auto& shape = basis.modes[m].coefficient[0];
                        state[6 * m + 3 * side + order] = inner(shape, field) / inner(shape, shape);
                    }
                }
            }
            const auto expected = fuelsim::evaluate_modal_beam(section, original, element, reference, false);
            const auto actual = fuelsim::evaluate_modal_beam(section, basis, element, state, false);
            require(std::abs(actual.energy / expected.energy - 1.0) < 1e-10,
                "Width enrichment changed classical extension, bending or torsion energy");
            for (std::size_t q = 0; q < actual.strain.size(); ++q)
                for (std::size_t c = 0; c < 6; ++c)
                    require(std::abs(actual.strain[q][c] - expected.strain[q][c]) < 1e-12,
                        "Width enrichment changed classical strain or shear cancellation");
        }
        if (strips == 2) {
            const auto size = 6 * basis.modes.size();
            const auto stiffness =
                fuelsim::modal_beam_linear_stiffness(fuelsim::linearize_modal_section(section, basis), element);
            std::vector<double> state(size), direction(size), plus(size), minus(size);
            constexpr double step = 1e-7;
            for (std::size_t i = 0; i < size; ++i) {
                state[i] = 1e-5 * std::sin(0.17 * static_cast<double>(i));
                direction[i] = std::cos(0.13 * static_cast<double>(i));
                plus[i] = state[i] + step * direction[i];
                minus[i] = state[i] - step * direction[i];
            }
            const auto actual = fuelsim::evaluate_modal_beam(section, basis, element, state, false);
            const auto rp = fuelsim::evaluate_modal_beam(section, basis, element, plus, false);
            const auto rm = fuelsim::evaluate_modal_beam(section, basis, element, minus, false);
            double error = 0.0, norm = 0.0, work = 0.0;
            for (std::size_t i = 0; i < size; ++i) {
                double exact = 0.0;
                for (std::size_t j = 0; j < size; ++j)
                    exact += stiffness[size * i + j] * direction[j];
                error = std::hypot(error, (rp.residual[i] - rm.residual[i]) / (2 * step) - exact);
                norm = std::hypot(norm, exact);
                work += state[i] * actual.residual[i];
            }
            std::cout << "width_directional_derivative_error=" << error / norm << '\n';
            require(error / norm < 1e-8 && std::abs(work / (2 * actual.energy) - 1.0) < 1e-10,
                "Width coupling failed the material-point tangent or virtual work");
        }
        previous = basis;
    }
    const auto asymmetric = fuelsim::enrich_section_width(section, seed, {-0.01, -0.0025, 0.01});
    require(!asymmetric.reflected_nodes[0].empty() && asymmetric.reflected_nodes[1].empty(),
        "An asymmetric width partition acquired unrequested reflected line fields");
    bool rejected = false;
    try {
        (void)fuelsim::enrich_section_width(section, seed, {-0.01, -0.003, 0.01});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "A width boundary cutting a section element was accepted");
}

void verify_condensation() {
    const auto section =
        fuelsim::test::rectangle({-0.001, 0.001}, 0.02, {fuelsim::test::region("al", 70e9, 0.3)}, false, 2);
    const auto basis = fuelsim::build_reduced_section_basis(section, false, 0);
    const auto tangent = fuelsim::linearize_modal_section(section, basis);
    // Four physical quintic beam elements, both end regions eliminated. Root
    // jets use multipliers so the eliminated block includes an indefinite
    // constraint system and a nonzero prescribed-displacement right hand side.
    constexpr PetscInt nodes = 5, amplitudes = 45, constraints = 9, total = amplitudes + constraints;
    Mat matrix = nullptr;
    Vec exact = nullptr, rhs = nullptr, solution = nullptr, work = nullptr;
    const auto check = [](PetscErrorCode code) {
        require(code == PETSC_SUCCESS, "Condensation test PETSc failure");
    };
    check(MatCreateAIJ(PETSC_COMM_WORLD,
        PETSC_DECIDE,
        PETSC_DECIDE,
        total,
        total,
        total,
        nullptr,
        total,
        nullptr,
        &matrix));
    check(MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    if (PetscGlobalRank == 0) {
        for (PetscInt e = 0; e < nodes - 1; ++e) {
            const auto stiffness = fuelsim::modal_beam_linear_stiffness(tangent,
                fuelsim::make_modal_beam_element(0.05 * static_cast<double>(e), 0.05 * static_cast<double>(e + 1)));
            std::array<PetscInt, 18> ids{};
            for (PetscInt m = 0; m < 3; ++m)
                for (PetscInt j = 0; j < 6; ++j)
                    ids[static_cast<std::size_t>(6 * m + j)] = (3 * m + j % 3) * nodes + e + j / 3;
            check(MatSetValues(matrix, 18, ids.data(), 18, ids.data(), stiffness.data(), ADD_VALUES));
        }
        for (PetscInt c = 0; c < constraints; ++c) {
            check(MatSetValue(matrix, amplitudes + c, c * nodes, 1.0, ADD_VALUES));
            check(MatSetValue(matrix, c * nodes, amplitudes + c, 1.0, ADD_VALUES));
        }
    }
    check(MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY));
    check(MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY));
    check(MatCreateVecs(matrix, &exact, &rhs));
    check(VecDuplicate(exact, &solution));
    check(VecDuplicate(exact, &work));
    // Equilibrate physical stiffness and multiplier rows without changing the
    // underlying virtual-work identity or removing any coupling.
    check(MatGetDiagonal(matrix, work));
    PetscInt begin, end;
    check(VecGetOwnershipRange(work, &begin, &end));
    PetscScalar* scale;
    check(VecGetArray(work, &scale));
    for (PetscInt i = begin; i < end; ++i)
        scale[i - begin] = i < amplitudes ? 1.0 / std::sqrt(scale[i - begin]) : 1.0;
    check(VecRestoreArray(work, &scale));
    check(MatDiagonalScale(matrix, work, work));
    std::vector<bool> eliminated(static_cast<std::size_t>(total), true);
    for (PetscInt i = 0; i < amplitudes; ++i)
        eliminated[static_cast<std::size_t>(i)] = i % nodes < 2 || i % nodes == 4;
    {
        fuelsim::ModalStaticCondensation condensed(matrix, eliminated);
        for (int load = 1; load <= 2; ++load) {
            check(VecGetArray(exact, &scale));
            for (PetscInt i = begin; i < end; ++i)
                scale[i - begin] = std::sin(0.17 * load * static_cast<double>(i + 1));
            check(VecRestoreArray(exact, &scale));
            check(MatMult(matrix, exact, rhs));
            condensed.solve(rhs, solution);
            PetscScalar expected_work, actual_work;
            check(VecDot(exact, rhs, &expected_work));
            check(VecDot(solution, rhs, &actual_work));
            require(std::abs(actual_work / expected_work - 1.0) < 1e-9, "Condensation changed physical virtual work");
            check(VecWAXPY(work, -1.0, exact, solution));
            PetscReal error, norm;
            check(VecNorm(work, NORM_2, &error));
            check(VecNorm(exact, NORM_2, &norm));
            std::cout << "condensation_solution_error=" << error / norm << '\n';
            require(error / norm < 1e-8, "Condensed and recovered beam solution differs from the original system");
            check(MatMult(matrix, solution, work));
            check(VecAXPY(work, -1.0, rhs));
            check(VecNorm(work, NORM_2, &error));
            check(VecNorm(rhs, NORM_2, &norm));
            require(error / norm < 1e-10, "Condensed beam failed the original equations");
        }
    }
    check(VecDestroy(&work));
    check(VecDestroy(&solution));
    check(VecDestroy(&rhs));
    check(VecDestroy(&exact));
    check(MatDestroy(&matrix));
}

void verify() {
    for (std::size_t side = 0; side < 2; ++side) {
        const double lower = 0.027, upper = 0.027312729;
        const auto shape = fuelsim::modal_beam_shape(lower, upper, side == 0 ? lower : upper);
        for (std::size_t order = 0; order < 3; ++order)
            for (std::size_t column = 0; column < 6; ++column)
                require(shape[order][column] == (column == 3 * side + order ? 1.0 : 0.0),
                    "Short-element endpoint values must be the exact nodal Hermite jets");
    }
    const auto section =
        fuelsim::test::rectangle({-0.001, 0.0, 0.001}, 0.02, {fuelsim::test::region("al", 70e9, 0.3)}, false, 2);
    const auto basis = fuelsim::build_reduced_section_basis(section, true, 24);
    const auto extended = fuelsim::build_reduced_section_basis(section, true, 28);
    // Truncation must preserve each actual mirror symmetry. Otherwise a
    // homogeneous axial load can acquire spurious bending even when every
    // retained field is independently normalized and the solve is accurate.
    for (std::size_t axis = 0; axis < 2; ++axis) {
        std::vector<std::size_t> reflected;
        for (const auto& node : section.nodes()) {
            const auto found = std::find_if(section.nodes().begin(), section.nodes().end(), [&](const auto& other) {
                return std::abs(other.x - (axis == 0 ? -node.x : node.x)) < 1e-14
                       && std::abs(other.y - (axis == 1 ? -node.y : node.y)) < 1e-14;
            });
            require(found != section.nodes().end(), "Symmetric section fixture is missing a reflected node");
            reflected.push_back(static_cast<std::size_t>(found - section.nodes().begin()));
        }
        for (const auto& mode : extended.modes) {
            double even = 0.0, odd = 0.0, norm = 0.0;
            const auto n = section.nodes().size();
            for (std::size_t c = 0; c < 3; ++c)
                for (std::size_t i = 0; i < n; ++i) {
                    const double value = mode.coefficient[0][c * n + i];
                    const double mirror = (c == axis ? -1.0 : 1.0) * mode.coefficient[0][c * n + reflected[i]];
                    norm = std::hypot(norm, value);
                    even = std::hypot(even, value - mirror);
                    odd = std::hypot(odd, value + mirror);
                }
            require(norm > 0.0 && std::min(even, odd) / norm < 1e-9,
                "A truncated section direction mixes distinct mirror symmetry spaces");
        }
    }
    for (std::size_t i = 0; i < basis.modes.size(); ++i)
        require(basis.modes[i].coefficient == extended.modes[i].coefficient,
            "Increasing the requested mode count changed an already retained section direction");
    require(basis.modes[6].kind == fuelsim::SectionMode::Kind::poisson_relaxation
                && basis.modes[9].kind == fuelsim::SectionMode::Kind::axial_warping,
        "Expected independent Poisson and axial warping relaxation");
    const auto axial_inner = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double value = 0.0;
        for (const auto& point : section.points()) {
            double left = 0.0, right = 0.0;
            for (std::size_t i = 0; i < 8; ++i) {
                left += point.shape[i] * a[point.nodes[i]];
                right += point.shape[i] * b[point.nodes[i]];
            }
            value += point.weight * left * right;
        }
        return value;
    };
    std::vector<std::vector<double>> axial_fields;
    for (std::size_t global = 0; global < 3; ++global) {
        std::vector<double> field(section.nodes().size(), 1.0);
        for (std::size_t i = 0; i < field.size(); ++i) {
            if (global == 1)
                field[i] = section.nodes()[i].x;
            if (global == 2)
                field[i] = section.nodes()[i].y;
        }
        const double norm = std::sqrt(axial_inner(field, field));
        for (double& value : field)
            value /= norm;
        axial_fields.push_back(std::move(field));
    }
    for (const auto& mode : basis.modes)
        if (mode.kind == fuelsim::SectionMode::Kind::axial_warping) {
            const auto& coefficient = mode.coefficient[0];
            std::vector<double> field(coefficient.begin() + static_cast<std::ptrdiff_t>(2 * section.nodes().size()),
                coefficient.end());
            require(std::abs(axial_inner(field, field) - 1.0) < 1e-10, "Axial warping normalization failed");
            for (const auto& previous : axial_fields)
                require(std::abs(axial_inner(previous, field)) < 1e-10,
                    "Independent axial warping overlaps the classical or preceding axial space");
            axial_fields.push_back(std::move(field));
            for (std::size_t q = 0; q < section.points().size(); ++q) {
                const auto k = fuelsim::modal_section_kinematics(section, mode, q, {1.0, 0.0, 0.0, 0.0});
                require(k.displacement[0] == 0.0 && k.displacement[1] == 0.0 && k.gradient[8] == 0.0,
                    "Independent constant-amplitude axial warping has incorrect kinematics");
            }
        }
    // A physical fixed end must admit finite axial strain and bending curvature
    // once independent Poisson relaxation fields are present.
    for (std::size_t mode = 0; mode < 3; ++mode) {
        const auto order = mode == 0 ? 1U : 2U;
        double product = 0.0, norm = 0.0;
        for (std::size_t i = 0; i < 2 * section.nodes().size(); ++i) {
            const double correction = basis.modes[mode].coefficient[order][i];
            const double relaxation = basis.modes[6 + mode].coefficient[0][i];
            product += correction * relaxation;
            norm += relaxation * relaxation;
        }
        std::array<double, 4> amplitude{};
        amplitude[order] = 1.0;
        for (std::size_t q = 0; q < section.points().size(); ++q) {
            const auto classic = fuelsim::modal_section_kinematics(section, basis.modes[mode], q, amplitude);
            const auto relaxation =
                fuelsim::modal_section_kinematics(section, basis.modes[6 + mode], q, {-product / norm, 0.0, 0.0, 0.0});
            for (std::size_t c = 0; c < 3; ++c)
                require(std::abs(classic.displacement[c] + relaxation.displacement[c]) < 1.0e-15,
                    "Physical clamp spuriously suppresses strain or curvature");
        }
    }
    const auto element = fuelsim::make_modal_beam_element(0.0, 0.5);
    const auto n = 6 * basis.modes.size();
    std::vector<double> state(n);
    for (std::size_t direction = 0; direction < 2; ++direction) {
        state.assign(n, 0.0);
        state[6 * (4 + direction)] = state[6 * (4 + direction) + 3] = 0.001;
        const auto response = fuelsim::evaluate_modal_beam(section, basis, element, state);
        const double expected = 0.5 * (70e9 / 2.6) * section.area() * 0.5 * 1e-6;
        require(std::abs(response.energy / expected - 1.0) < 1e-10,
            "Independent rotation did not reproduce transverse shear energy");
        for (const auto& strain : response.strain)
            require(std::abs(strain[direction == 0 ? 5 : 4] - 0.0005) < 1e-15,
                "Independent rotation shear kinematics failed");
    }
    state.assign(n, 0.0);
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
    const auto mixed_basis = fuelsim::make_mixed_bending_basis(basis);
    const auto independent_basis = fuelsim::make_independent_section_basis(section, mixed_basis);
    auto mixed_state = state;
    mixed_state[6 * 4 + 1] = mixed_state[6 * 4 + 4] = -0.02;
    mixed_state[6 * 4 + 3] = -0.02 * 0.5;
    const auto mixed_bending = fuelsim::evaluate_modal_beam(section, mixed_basis, element, mixed_state);
    require(std::abs(mixed_bending.energy / bending.energy - 1.0) < 1e-10,
        "Independent rotations lost the classical pure bending energy");
    for (std::size_t q = 0; q < bending.strain.size(); ++q)
        for (std::size_t c = 0; c < 6; ++c)
            require(std::abs(mixed_bending.strain[q][c] - bending.strain[q][c]) < 1e-15,
                "Independent rotations lost the Euler-Bernoulli strain or shear cancellation");
    // In this symmetric fixture Phi0 fields are mutually L2 orthogonal. Project
    // four known classical displacement polynomials and their derivatives onto
    // the independent physical amplitudes, then compare the full point strains.
    const auto displacement_inner = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double value = 0.0;
        const auto nodes = section.nodes().size();
        for (const auto& point : section.points())
            for (std::size_t c = 0; c < 3; ++c) {
                double left = 0.0, right = 0.0;
                for (std::size_t i = 0; i < 8; ++i) {
                    left += point.shape[i] * a[c * nodes + point.nodes[i]];
                    right += point.shape[i] * b[c * nodes + point.nodes[i]];
                }
                value += point.weight * left * right;
            }
        return value;
    };
    for (std::size_t classic = 0; classic < 4; ++classic) {
        std::vector<double> original(n), physical(n);
        for (std::size_t side = 0; side < 2; ++side) {
            const double z = side == 0 ? element.lower : element.upper;
            const std::array<double, 5> derivative =
                classic == 0 || classic == 3 ? std::array<double, 5>{0.001 * z, 0.001, 0.0, 0.0, 0.0}
                                             : std::array<double, 5>{0.01 * z * z, 0.02 * z, 0.02, 0.0, 0.0};
            for (std::size_t order = 0; order < 3; ++order) {
                original[6 * classic + 3 * side + order] = derivative[order];
                std::vector<double> field(3 * section.nodes().size());
                for (std::size_t k = 0; k < 3; ++k)
                    for (std::size_t i = 0; i < field.size(); ++i)
                        field[i] += basis.modes[classic].coefficient[k][i] * derivative[order + k];
                for (std::size_t mode = 0; mode < independent_basis.modes.size(); ++mode) {
                    const auto& shape = independent_basis.modes[mode].coefficient[0];
                    physical[6 * mode + 3 * side + order] =
                        displacement_inner(shape, field) / displacement_inner(shape, shape);
                }
            }
        }
        const auto expected = fuelsim::evaluate_modal_beam(section, basis, element, original, false);
        const auto actual = fuelsim::evaluate_modal_beam(section, independent_basis, element, physical, false);
        auto minimum_basis = independent_basis;
        minimum_basis.modes.resize(10);
        const std::vector<double> minimum_state(physical.begin(), physical.begin() + 60);
        const auto minimum = fuelsim::evaluate_modal_beam(section, minimum_basis, element, minimum_state, false);
        require(std::abs(minimum.energy / expected.energy - 1.0) < 1e-10,
            "The minimum retained section basis lost a classical mode energy");
        require(std::abs(actual.energy / expected.energy - 1.0) < 1e-10,
            "Independent physical amplitudes lost a classical extension, bending or torsion energy");
        for (std::size_t q = 0; q < actual.strain.size(); ++q)
            for (std::size_t c = 0; c < 6; ++c)
                require(std::abs(actual.strain[q][c] - expected.strain[q][c]) < 1e-13,
                    "Independent physical amplitudes changed a classical strain or shear cancellation");
        for (std::size_t q = 0; q < minimum.strain.size(); ++q)
            for (std::size_t c = 0; c < 6; ++c)
                require(std::abs(minimum.strain[q][c] - expected.strain[q][c]) < 1e-13,
                    "The minimum retained section basis lost a classical point strain");
    }
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
        if (rigid == 3 || rigid == 4) {
            const auto rotation = 6 * (rigid + 1);
            state[rotation] = state[rotation + 3] = -1.0;
        }
        const auto mixed = fuelsim::evaluate_modal_beam(section, mixed_basis, element, state, false);
        require(mixed.energy < 1e-16, "Mixed beam rigid motion acquired strain energy");
    }
    for (const auto* evaluated_basis : {&basis, &mixed_basis, &independent_basis}) {
        for (std::size_t i = 0; i < n; ++i)
            state[i] = 1.0e-4 * std::sin(static_cast<double>(i + 1));
        const auto response = fuelsim::evaluate_modal_beam(section, *evaluated_basis, element, state);
        const auto sampled = fuelsim::sample_modal_section(section, *evaluated_basis);
        const auto cached = fuelsim::evaluate_modal_beam(section, *evaluated_basis, element, state, true, &sampled);
        require(cached.residual == response.residual && cached.jacobian == response.jacobian
                    && cached.strain == response.strain && cached.stress == response.stress
                    && cached.energy == response.energy,
            "Reusing fixed section kinematics changed the material-point response or tangent");
        const auto integrated =
            fuelsim::modal_beam_linear_stiffness(fuelsim::linearize_modal_section(section, *evaluated_basis), element);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                require(std::abs(integrated[n * i + j] - response.jacobian[n * i + j])
                            < 1e-10 * std::sqrt(response.jacobian[n * i + i] * response.jacobian[n * j + j]),
                    "Separated linear section/axial integration changed the material-point tangent");
        const auto residual_only = fuelsim::evaluate_modal_beam(section, *evaluated_basis, element, state, false);
        require(response.residual == residual_only.residual, "Beam residual-only and tangent paths disagree");
        std::vector<double> direction(n), plus = state, minus = state;
        constexpr double step = 1.0e-7;
        for (std::size_t i = 0; i < n; ++i) {
            direction[i] = std::cos(static_cast<double>(i + 1));
            plus[i] += step * direction[i];
            minus[i] -= step * direction[i];
        }
        const auto rp = fuelsim::evaluate_modal_beam(section, *evaluated_basis, element, plus, false);
        const auto rm = fuelsim::evaluate_modal_beam(section, *evaluated_basis, element, minus, false);
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
        std::cout << "beam_mixed=" << (evaluated_basis == &mixed_basis)
                  << " directional_derivative_error=" << error / norm << " symmetry=" << symmetry << '\n';
        require(error / norm < 1.0e-8 && symmetry < 1.0e-10, "Beam consistent tangent failed");
        require(std::abs(work / response.energy - 1.0) < 1.0e-10, "Beam energy and residual virtual work disagree");
        require(response.strain.size() == 6 * section.points().size(), "Section material points were reduced away");
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "Modal beam local contracts");
        verify_reflection_tangent();
        verify_solid_interface();
        verify_width_enrichment();
        verify_condensation();
        verify();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
