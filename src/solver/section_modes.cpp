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

class SectionRelaxationSolver::Impl final {
  public:
    explicit Impl(CrossSection section) : _section(std::move(section)) {
        PetscBool initialized = PETSC_FALSE;
        check(PetscInitialized(&initialized));
        if (!initialized)
            throw std::logic_error("Section preprocessing requires a PETSc session");
        const auto n = _section.nodes().size();
        if (n > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()) / 2)
            throw std::invalid_argument("Section exceeds PETSc index capacity");
        const auto size = static_cast<PetscInt>(2 * n);
        check(MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 32, nullptr, &_system.matrix));
        check(MatSetOption(_system.matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
        check(VecCreateSeq(PETSC_COMM_SELF, size, &_system.rhs));
        check(VecDuplicate(_system.rhs, &_system.solution));
        check(VecDuplicate(_system.rhs, &_system.defect));
        constexpr std::array<double, 6> metric{1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
        for (const auto& point : _section.points()) {
            const auto& region = _section.regions()[point.region];
            const auto material = evaluate_stress_tangent(region.material,
                {},
                region.temperature,
                0.0,
                nullptr,
                material_context(0.0, point.position));
            std::array<double, 6> coupling{};
            for (std::size_t a = 0; a < 6; ++a)
                coupling[a] = metric[a] * material.tangent[a][2];
            _coupling.push_back(coupling);
            _shear.push_back(0.5 * material.tangent[3][3]);
            std::array<PetscInt, 16> ids{};
            std::array<PetscScalar, 256> stiffness{};
            for (std::size_t i = 0; i < 16; ++i) {
                ids[i] = static_cast<PetscInt>(point.nodes[i % 8] + (i / 8) * n);
                const auto bi = transverse_column(point, i);
                for (std::size_t a = 0; a < 6; ++a)
                    for (std::size_t j = 0; j < 16; ++j) {
                        const auto bj = transverse_column(point, j);
                        for (std::size_t b = 0; b < 6; ++b)
                            stiffness[16 * i + j] += point.weight * metric[a] * bi[a] * material.tangent[a][b] * bj[b];
                    }
            }
            check(MatSetValues(_system.matrix, 16, ids.data(), 16, ids.data(), stiffness.data(), ADD_VALUES));
        }
        check(MatAssemblyBegin(_system.matrix, MAT_FINAL_ASSEMBLY));
        check(MatAssemblyEnd(_system.matrix, MAT_FINAL_ASSEMBLY));
        check(MatDuplicate(_system.matrix, MAT_COPY_VALUES, &_system.unconstrained));
        // Remove exactly two translations and one rotation, then restore the
        // integral L2 gauge after solving and check the original equation.
        std::size_t far = 0;
        double separation = 0.0;
        bool constrain_y = true;
        for (std::size_t i = 1; i < n; ++i) {
            const double dx = std::abs(_section.nodes()[i].x - _section.nodes()[0].x);
            const double dy = std::abs(_section.nodes()[i].y - _section.nodes()[0].y);
            if (std::max(dx, dy) > separation) {
                separation = std::max(dx, dy);
                far = i;
                constrain_y = dx >= dy;
            }
        }
        _pins = {0, static_cast<PetscInt>(n), static_cast<PetscInt>(far + (constrain_y ? n : 0))};
        check(MatZeroRowsColumns(_system.matrix, 3, _pins.data(), 1.0, nullptr, nullptr));
        check(KSPCreate(PETSC_COMM_SELF, &_system.solver));
        check(KSPSetOperators(_system.solver, _system.matrix, _system.matrix));
        check(KSPSetType(_system.solver, KSPPREONLY));
        PC pc = nullptr;
        check(KSPGetPC(_system.solver, &pc));
        check(PCSetType(pc, PCLU));
        check(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    }

    SectionTransverseRelaxation solve(const std::vector<double>& axial_strain) const {
        const auto n = _section.nodes().size();
        if (axial_strain.size() != n)
            throw std::invalid_argument("Section relaxation requires one axial strain per section node");
        for (double value : axial_strain)
            if (!std::isfinite(value))
                throw std::invalid_argument("Section relaxation requires finite axial strain");
        std::vector<double> load(2 * n);
        for (std::size_t q = 0; q < _section.points().size(); ++q) {
            const auto& point = _section.points()[q];
            double axial = axial_strain[point.nodes[0]];
            for (std::size_t i = 0; i < 8; ++i)
                axial += point.shape[i] * (axial_strain[point.nodes[i]] - axial_strain[point.nodes[0]]);
            for (std::size_t i = 0; i < 16; ++i) {
                const auto bi = transverse_column(point, i);
                for (std::size_t a = 0; a < 6; ++a)
                    load[point.nodes[i % 8] + (i / 8) * n] -= point.weight * bi[a] * _coupling[q][a] * axial;
            }
        }
        return solve_load(std::move(load));
    }

    SectionTransverseRelaxation solve_transverse_corrector(const std::vector<double>& source) const {
        const auto n = _section.nodes().size();
        if (source.size() != 2 * n)
            throw std::invalid_argument("Transverse corrector requires field-major [ux,uy] values");
        for (double value : source)
            if (!std::isfinite(value))
                throw std::invalid_argument("Transverse corrector requires finite values");
        std::vector<double> load(2 * n);
        for (std::size_t q = 0; q < _section.points().size(); ++q) {
            const auto& point = _section.points()[q];
            for (std::size_t c = 0; c < 2; ++c) {
                double value = 0.0;
                for (std::size_t i = 0; i < 8; ++i)
                    value += point.shape[i] * source[c * n + point.nodes[i]];
                for (std::size_t i = 0; i < 8; ++i)
                    load[c * n + point.nodes[i]] += point.weight * _shear[q] * point.shape[i] * value;
            }
        }
        // For u_perp=p*q, 3D shear contributes (q')^2*integral(G*p.p)/2.
        // Transverse equilibrium supplies KD*chi=H*p. Global force/moment
        // belongs to the explicit beam fields: project this load into their
        // W-orthogonal complement before the three scalar gauges are pinned.
        const auto center = _section.centroid();
        std::array<double, 3> force{}, norm{};
        std::array<std::vector<double>, 3> weights;
        for (auto& weight : weights)
            weight.resize(2 * n);
        for (std::size_t i = 0; i < n; ++i) {
            force[0] += load[i];
            force[1] += load[n + i];
            force[2] +=
                -(_section.nodes()[i].y - center.y) * load[i] + (_section.nodes()[i].x - center.x) * load[n + i];
        }
        for (const auto& point : _section.points()) {
            const double x = point.position.x - center.x, y = point.position.y - center.y;
            norm[0] += point.weight;
            norm[1] += point.weight;
            norm[2] += point.weight * (x * x + y * y);
            for (std::size_t i = 0; i < 8; ++i) {
                const double weight = point.weight * point.shape[i];
                weights[0][point.nodes[i]] += weight;
                weights[1][n + point.nodes[i]] += weight;
                weights[2][point.nodes[i]] -= weight * y;
                weights[2][n + point.nodes[i]] += weight * x;
            }
        }
        for (std::size_t i = 0; i < 2 * n; ++i)
            for (std::size_t global = 0; global < 3; ++global)
                load[i] -= weights[global][i] * force[global] / norm[global];
        return solve_load(std::move(load));
    }

  private:
    SectionTransverseRelaxation solve_load(std::vector<double> load) const {
        const auto n = _section.nodes().size();
        PetscScalar* values = nullptr;
        check(VecGetArray(_system.rhs, &values));
        std::copy(load.begin(), load.end(), values);
        for (auto pin : _pins)
            values[pin] = 0.0;
        check(VecRestoreArray(_system.rhs, &values));
        check(KSPSolve(_system.solver, _system.rhs, _system.solution));
        KSPConvergedReason reason;
        check(KSPGetConvergedReason(_system.solver, &reason));
        if (reason <= 0)
            throw std::runtime_error("Section auxiliary solve did not converge");
        SectionTransverseRelaxation result;
        const PetscScalar* solution = nullptr;
        check(VecGetArrayRead(_system.solution, &solution));
        result.transverse.assign(solution, solution + static_cast<PetscInt>(2 * n));
        check(VecRestoreArrayRead(_system.solution, &solution));
        remove_rigid_motion(_section, result.transverse);
        check(VecGetArray(_system.solution, &values));
        std::copy(result.transverse.begin(), result.transverse.end(), values);
        check(VecRestoreArray(_system.solution, &values));
        check(VecGetArray(_system.rhs, &values));
        std::copy(load.begin(), load.end(), values);
        check(VecRestoreArray(_system.rhs, &values));
        check(MatMult(_system.unconstrained, _system.solution, _system.defect));
        check(VecAXPY(_system.defect, -1.0, _system.rhs));
        PetscReal error = 0.0, norm = 0.0;
        check(VecNorm(_system.defect, NORM_2, &error));
        check(VecNorm(_system.rhs, NORM_2, &norm));
        result.relative_residual = norm > 0.0 ? error / norm : error;
        if (!std::isfinite(result.relative_residual) || result.relative_residual > 1.0e-9)
            throw std::runtime_error("Section unconstrained equilibrium residual exceeds tolerance");
        return result;
    }

  private:
    CrossSection _section;
    AuxiliarySystem _system;
    std::array<PetscInt, 3> _pins{};
    std::vector<std::array<double, 6>> _coupling;
    std::vector<double> _shear;
};

SectionRelaxationSolver::SectionRelaxationSolver(CrossSection section)
    : _impl(std::make_unique<Impl>(std::move(section))) {
}

SectionRelaxationSolver::~SectionRelaxationSolver() = default;

SectionTransverseRelaxation SectionRelaxationSolver::solve(const std::vector<double>& axial_strain) const {
    return _impl->solve(axial_strain);
}

SectionTransverseRelaxation SectionRelaxationSolver::solve_transverse_corrector(
    const std::vector<double>& source) const {
    return _impl->solve_transverse_corrector(source);
}

ClassicSectionBasis build_classic_section_basis(const CrossSection& section) {
    SectionRelaxationSolver solver(section);
    ClassicSectionBasis basis{};
    basis.origin = section.elastic_center();
    for (std::size_t mode = 0; mode < 3; ++mode) {
        std::vector<double> axial(section.nodes().size(), 1.0);
        for (std::size_t i = 0; i < axial.size(); ++i) {
            if (mode == 1)
                axial[i] = -(section.nodes()[i].x - basis.origin.x);
            if (mode == 2)
                axial[i] = -(section.nodes()[i].y - basis.origin.y);
        }
        auto relaxation = solver.solve(axial);
        basis.modes[mode].kind = static_cast<ClassicSectionMode::Kind>(mode);
        basis.modes[mode].correction = std::move(relaxation.transverse);
        basis.modes[mode].auxiliary_relative_residual = relaxation.relative_residual;
    }
    return basis;
}
} // namespace fuelsim
