#include "solver/section_modes.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <petscksp.h>
#include <stdexcept>
#include <string>

namespace fuelsim {
namespace {
void check(PetscErrorCode code) {
    if (code != PETSC_SUCCESS)
        throw std::runtime_error("Section auxiliary PETSc failure: " + std::to_string(code));
}

struct AuxiliarySystem final {
    Mat matrix = nullptr, unconstrained = nullptr;
    Vec rhs = nullptr, solution = nullptr, defect = nullptr;
    KSP solver = nullptr;

    ~AuxiliarySystem() {
        (void)KSPDestroy(&solver);
        (void)VecDestroy(&defect);
        (void)VecDestroy(&solution);
        (void)VecDestroy(&rhs);
        (void)MatDestroy(&matrix);
        (void)MatDestroy(&unconstrained);
    }
};

// Tensor shear components: engineering shear factors enter the virtual-work metric.
std::array<double, 6> transverse_column(const SectionPoint& point, std::size_t local) {
    std::array<double, 6> b{};
    const auto& gradient = point.gradient[local % 8];
    if (local < 8) {
        b[0] = gradient[0];
        b[3] = 0.5 * gradient[1];
    } else {
        b[1] = gradient[1];
        b[3] = 0.5 * gradient[0];
    }
    return b;
}

void remove_rigid_motion(const CrossSection& section, std::vector<double>& field) {
    const auto center = section.centroid();
    const auto n = section.nodes().size();
    double mean_x = 0.0, mean_y = 0.0, rotation = 0.0, inertia = 0.0;
    for (const auto& point : section.points()) {
        double ux = 0.0, uy = 0.0;
        for (std::size_t i = 0; i < 8; ++i) {
            ux += point.shape[i] * field[point.nodes[i]];
            uy += point.shape[i] * field[n + point.nodes[i]];
        }
        const double x = point.position.x - center.x, y = point.position.y - center.y;
        mean_x += point.weight * ux;
        mean_y += point.weight * uy;
        rotation += point.weight * (-y * ux + x * uy);
        inertia += point.weight * (x * x + y * y);
    }
    mean_x /= section.area();
    mean_y /= section.area();
    rotation /= inertia;
    for (std::size_t i = 0; i < n; ++i) {
        field[i] -= mean_x - rotation * (section.nodes()[i].y - center.y);
        field[n + i] -= mean_y + rotation * (section.nodes()[i].x - center.x);
    }
}
} // namespace

ClassicSectionBasis build_classic_section_basis(const CrossSection& section) {
    PetscBool initialized = PETSC_FALSE;
    check(PetscInitialized(&initialized));
    if (!initialized)
        throw std::logic_error("Section preprocessing requires a PETSc session");
    const auto n = section.nodes().size();
    if (n > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()) / 2)
        throw std::invalid_argument("Section exceeds PETSc index capacity");
    const auto size = static_cast<PetscInt>(2 * n);
    AuxiliarySystem system;
    check(MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 32, nullptr, &system.matrix));
    check(MatSetOption(system.matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    check(VecCreateSeq(PETSC_COMM_SELF, size, &system.rhs));
    check(VecDuplicate(system.rhs, &system.solution));
    check(VecDuplicate(system.rhs, &system.defect));
    std::array<std::vector<double>, 3> loads;
    for (auto& load : loads)
        load.resize(2 * n);
    constexpr std::array<double, 6> metric{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    const auto origin = section.elastic_center();
    for (const auto& point : section.points()) {
        const auto& region = section.regions()[point.region];
        const auto material = evaluate_stress_tangent(region.material,
            {},
            region.temperature,
            0.0,
            nullptr,
            material_context(0.0, point.position));
        std::array<PetscInt, 16> ids{};
        std::array<PetscScalar, 256> stiffness{};
        const std::array<double, 3> axial{1.0, -(point.position.x - origin.x), -(point.position.y - origin.y)};
        for (std::size_t i = 0; i < 16; ++i) {
            ids[i] = static_cast<PetscInt>(point.nodes[i % 8] + (i / 8) * n);
            const auto bi = transverse_column(point, i);
            for (std::size_t a = 0; a < 6; ++a) {
                for (std::size_t mode = 0; mode < 3; ++mode)
                    loads[mode][static_cast<std::size_t>(ids[i])] -=
                        point.weight * metric[a] * bi[a] * material.tangent[a][2] * axial[mode];
                for (std::size_t j = 0; j < 16; ++j) {
                    const auto bj = transverse_column(point, j);
                    for (std::size_t b = 0; b < 6; ++b)
                        stiffness[16 * i + j] += point.weight * metric[a] * bi[a] * material.tangent[a][b] * bj[b];
                }
            }
        }
        check(MatSetValues(system.matrix, 16, ids.data(), 16, ids.data(), stiffness.data(), ADD_VALUES));
    }
    check(MatAssemblyBegin(system.matrix, MAT_FINAL_ASSEMBLY));
    check(MatAssemblyEnd(system.matrix, MAT_FINAL_ASSEMBLY));
    check(MatDuplicate(system.matrix, MAT_COPY_VALUES, &system.unconstrained));
    // Three independent gauge constraints remove exactly the two translations and one rotation.
    // Afterwards project to the integral L2-orthogonal gauge, important for axial derivatives.
    std::size_t far = 0;
    double separation = 0.0;
    bool constrain_y = true;
    for (std::size_t i = 1; i < n; ++i) {
        const double dx = std::abs(section.nodes()[i].x - section.nodes()[0].x);
        const double dy = std::abs(section.nodes()[i].y - section.nodes()[0].y);
        if (std::max(dx, dy) > separation) {
            separation = std::max(dx, dy);
            far = i;
            constrain_y = dx >= dy;
        }
    }
    const std::array<PetscInt, 3> pins{0, static_cast<PetscInt>(n), static_cast<PetscInt>(far + (constrain_y ? n : 0))};
    check(MatZeroRowsColumns(system.matrix, 3, pins.data(), 1.0, nullptr, nullptr));
    check(KSPCreate(PETSC_COMM_SELF, &system.solver));
    check(KSPSetOperators(system.solver, system.matrix, system.matrix));
    check(KSPSetType(system.solver, KSPPREONLY));
    PC pc = nullptr;
    check(KSPGetPC(system.solver, &pc));
    check(PCSetType(pc, PCLU));
    check(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    ClassicSectionBasis basis{};
    basis.origin = origin;
    for (std::size_t mode = 0; mode < 3; ++mode) {
        auto& result = basis.modes[mode];
        result.kind = static_cast<ClassicSectionMode::Kind>(mode);
        auto constrained_load = loads[mode];
        for (auto pin : pins)
            constrained_load[static_cast<std::size_t>(pin)] = 0.0;
        PetscScalar* values = nullptr;
        check(VecGetArray(system.rhs, &values));
        std::copy(constrained_load.begin(), constrained_load.end(), values);
        check(VecRestoreArray(system.rhs, &values));
        check(KSPSolve(system.solver, system.rhs, system.solution));
        KSPConvergedReason reason;
        check(KSPGetConvergedReason(system.solver, &reason));
        if (reason <= 0)
            throw std::runtime_error("Section auxiliary solve did not converge");
        check(MatMult(system.matrix, system.solution, system.defect));
        check(VecAXPY(system.defect, -1.0, system.rhs));
        PetscReal error = 0.0, norm = 0.0;
        check(VecNorm(system.defect, NORM_2, &error));
        check(VecNorm(system.rhs, NORM_2, &norm));
        result.auxiliary_relative_residual = norm > 0.0 ? error / norm : error;
        if (result.auxiliary_relative_residual > 1.0e-9)
            throw std::runtime_error("Section auxiliary equation residual exceeds tolerance");
        const PetscScalar* solution = nullptr;
        check(VecGetArrayRead(system.solution, &solution));
        result.correction.assign(solution, solution + size);
        check(VecRestoreArrayRead(system.solution, &solution));
        remove_rigid_motion(section, result.correction);
        // Verify the original, unconstrained variational equation after fixing the L2 gauge.
        check(VecGetArray(system.solution, &values));
        std::copy(result.correction.begin(), result.correction.end(), values);
        check(VecRestoreArray(system.solution, &values));
        check(VecGetArray(system.rhs, &values));
        std::copy(loads[mode].begin(), loads[mode].end(), values);
        check(VecRestoreArray(system.rhs, &values));
        check(MatMult(system.unconstrained, system.solution, system.defect));
        check(VecAXPY(system.defect, -1.0, system.rhs));
        check(VecNorm(system.defect, NORM_2, &error));
        check(VecNorm(system.rhs, NORM_2, &norm));
        result.auxiliary_relative_residual = norm > 0.0 ? error / norm : error;
        if (result.auxiliary_relative_residual > 1.0e-9)
            throw std::runtime_error("Section unconstrained equilibrium residual exceeds tolerance");
    }
    return basis;
}
} // namespace fuelsim
