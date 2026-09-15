#include "solver/section_warping.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <petscksp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace fuelsim {
namespace {
void checked(PetscErrorCode code) {
    if (code != PETSC_SUCCESS)
        throw std::runtime_error("Section warping PETSc failure: " + std::to_string(code));
}

struct WarpingObjects final {
    Mat matrix = nullptr;
    Vec rhs = nullptr, solution = nullptr;
    KSP solver = nullptr;

    ~WarpingObjects() {
        (void)KSPDestroy(&solver);
        (void)VecDestroy(&solution);
        (void)VecDestroy(&rhs);
        (void)MatDestroy(&matrix);
    }
};

double shear_modulus(const CrossSection& section, const SectionPoint& point) {
    const auto& region = section.regions()[point.region];
    return region.material.active_properties(region.temperature, material_context(0.0, point.position))
        .shear_modulus.value();
}

void enforce_neumann_compatibility(std::vector<double>& load) {
    // The constant test function has identically zero gradient. Enforce this
    // exact assembly identity before the scalar gauge is pinned. In particular,
    // a solenoidal transverse field can have an almost zero load; a roundoff
    // component in the constant nullspace must not become a fictitious pin force.
    for (int pass = 0; pass < 2; ++pass) {
        double total = 0.0;
        for (double value : load)
            total += value;
        const double mean = total / static_cast<double>(load.size());
        for (double& value : load)
            value -= mean;
    }
}
} // namespace

class SectionWarpingSolver::Impl final {
  public:
    explicit Impl(CrossSection section) : _section(std::move(section)) {
        PetscBool initialized = PETSC_FALSE;
        checked(PetscInitialized(&initialized));
        if (!initialized)
            throw std::logic_error("Section warping requires a PETSc session");
        const auto n = _section.nodes().size();
        if (n > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
            throw std::invalid_argument("Warping mesh exceeds PETSc index capacity");
        const auto size = static_cast<PetscInt>(n);
        checked(MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, 24, nullptr, &_objects.matrix));
        checked(MatSetOption(_objects.matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
        checked(VecCreateSeq(PETSC_COMM_SELF, size, &_objects.rhs));
        checked(VecDuplicate(_objects.rhs, &_objects.solution));
        for (const auto& point : _section.points()) {
            std::array<PetscInt, 8> ids{};
            std::array<PetscScalar, 64> stiffness{};
            const double weight = point.weight * shear_modulus(_section, point);
            for (std::size_t i = 0; i < 8; ++i) {
                ids[i] = static_cast<PetscInt>(point.nodes[i]);
                for (std::size_t j = 0; j < 8; ++j)
                    stiffness[8 * i + j] =
                        weight
                        * (point.gradient[i][0] * point.gradient[j][0] + point.gradient[i][1] * point.gradient[j][1]);
            }
            checked(MatSetValues(_objects.matrix, 8, ids.data(), 8, ids.data(), stiffness.data(), ADD_VALUES));
        }
        checked(MatAssemblyBegin(_objects.matrix, MAT_FINAL_ASSEMBLY));
        checked(MatAssemblyEnd(_objects.matrix, MAT_FINAL_ASSEMBLY));
        const PetscInt pin = 0;
        checked(MatZeroRowsColumns(_objects.matrix, 1, &pin, 1.0, nullptr, nullptr));
        checked(KSPCreate(PETSC_COMM_SELF, &_objects.solver));
        checked(KSPSetOperators(_objects.solver, _objects.matrix, _objects.matrix));
        checked(KSPSetType(_objects.solver, KSPPREONLY));
        PC pc = nullptr;
        checked(KSPGetPC(_objects.solver, &pc));
        checked(PCSetType(pc, PCLU));
        checked(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    }

    SectionWarping solve(const std::vector<double>& transverse) const {
        const auto n = _section.nodes().size();
        if (transverse.size() != 2 * n)
            throw std::invalid_argument("Warping requires field-major [ux,uy] transverse values");
        for (double value : transverse)
            if (!std::isfinite(value))
                throw std::invalid_argument("Warping transverse field must be finite");
        std::vector<double> load(n, 0.0);
        for (const auto& point : _section.points()) {
            double ux = 0.0, uy = 0.0;
            for (std::size_t i = 0; i < 8; ++i) {
                ux += point.shape[i] * transverse[point.nodes[i]];
                uy += point.shape[i] * transverse[n + point.nodes[i]];
            }
            const double weight = point.weight * shear_modulus(_section, point);
            for (std::size_t i = 0; i < 8; ++i)
                load[point.nodes[i]] -= weight * (point.gradient[i][0] * ux + point.gradient[i][1] * uy);
        }
        enforce_neumann_compatibility(load);
        SectionWarping result;
        result.axial = solve_load(load);
        result.condensed_force.resize(2 * n);
        auto residual = load;
        for (double& value : residual)
            value = -value;
        // Apply the Schur complement through its minimizing field, never an explicit inverse.
        for (const auto& point : _section.points()) {
            double ux = 0.0, uy = 0.0, wx = 0.0, wy = 0.0;
            for (std::size_t i = 0; i < 8; ++i) {
                ux += point.shape[i] * transverse[point.nodes[i]];
                uy += point.shape[i] * transverse[n + point.nodes[i]];
                wx += point.gradient[i][0] * (result.axial[point.nodes[i]] - result.axial[point.nodes[0]]);
                wy += point.gradient[i][1] * (result.axial[point.nodes[i]] - result.axial[point.nodes[0]]);
            }
            const double gx = ux + wx, gy = uy + wy;
            const double weight = point.weight * shear_modulus(_section, point);
            result.stiffness += weight * (gx * gx + gy * gy);
            for (std::size_t i = 0; i < 8; ++i) {
                // Check A*omega-load using the already assembled load. For a
                // divergence-free transverse field its load can be nearly zero;
                // reintegrating grad(N)*(phi+grad(omega)) would subtract large
                // phi terms a second time and report their roundoff as solve error.
                residual[point.nodes[i]] += weight * (point.gradient[i][0] * wx + point.gradient[i][1] * wy);
                result.condensed_force[point.nodes[i]] += weight * point.shape[i] * gx;
                result.condensed_force[n + point.nodes[i]] += weight * point.shape[i] * gy;
            }
        }
        double defect = 0.0, norm = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            defect = std::hypot(defect, residual[i]);
            norm = std::hypot(norm, load[i]);
        }
        result.relative_residual = norm > 0.0 ? defect / norm : defect;
        if (!std::isfinite(result.relative_residual) || result.relative_residual > 1.0e-9) {
            std::ostringstream message;
            message << "Warping original equation residual exceeds tolerance: relative=" << result.relative_residual
                    << ", defect=" << defect << ", load=" << norm;
            throw std::runtime_error(message.str());
        }
        return result;
    }

    SectionAxialCorrector solve_axial_corrector(const std::vector<double>& source) const {
        const auto n = _section.nodes().size();
        if (source.size() != n)
            throw std::invalid_argument("Axial corrector requires one scalar source per section node");
        for (double value : source)
            if (!std::isfinite(value))
                throw std::invalid_argument("Axial corrector source must be finite");
        // For uz=psi*r, the pure axial terms of the 3D energy are
        // (r^2*psi^T*A_omega*psi + r'^2*psi^T*M_zz*psi)/2.
        // Axial equilibrium therefore supplies the next section direction
        // A_omega*chi=M_zz*psi. Subtract the Czz-weighted constant to satisfy
        // the natural Neumann compatibility condition. This generates higher
        // axial shapes from physical seeds without prescribing a polynomial.
        std::vector<double> value(_section.points().size()), weight(value.size());
        double mean = 0.0, total_weight = 0.0;
        for (std::size_t q = 0; q < value.size(); ++q) {
            const auto& point = _section.points()[q];
            const auto& region = _section.regions()[point.region];
            const auto properties =
                region.material.active_properties(region.temperature, material_context(0.0, point.position));
            weight[q] = point.weight * (properties.lame_lambda.value() + 2.0 * properties.shear_modulus.value());
            for (std::size_t i = 0; i < 8; ++i)
                value[q] += point.shape[i] * (source[point.nodes[i]] - source[0]);
            mean += weight[q] * value[q];
            total_weight += weight[q];
        }
        mean /= total_weight;
        std::vector<double> load(n);
        for (std::size_t q = 0; q < value.size(); ++q) {
            const auto& point = _section.points()[q];
            for (std::size_t i = 0; i < 8; ++i)
                load[point.nodes[i]] += weight[q] * point.shape[i] * (value[q] - mean);
        }
        enforce_neumann_compatibility(load);
        SectionAxialCorrector result;
        result.axial = solve_load(load);
        auto residual = load;
        for (double& entry : residual)
            entry = -entry;
        for (const auto& point : _section.points()) {
            std::array<double, 2> gradient{};
            for (std::size_t i = 0; i < 8; ++i)
                for (std::size_t c = 0; c < 2; ++c)
                    gradient[c] += point.gradient[i][c] * result.axial[point.nodes[i]];
            const double measure = point.weight * shear_modulus(_section, point);
            for (std::size_t i = 0; i < 8; ++i)
                residual[point.nodes[i]] +=
                    measure * (point.gradient[i][0] * gradient[0] + point.gradient[i][1] * gradient[1]);
        }
        double defect = 0.0, norm = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            defect = std::hypot(defect, residual[i]);
            norm = std::hypot(norm, load[i]);
        }
        result.relative_residual = norm > 0.0 ? defect / norm : defect;
        if (!std::isfinite(result.relative_residual) || result.relative_residual > 1e-9)
            throw std::runtime_error("Axial corrector original equation residual exceeds tolerance");
        return result;
    }

  private:
    std::vector<double> solve_load(const std::vector<double>& load) const {
        PetscScalar* values = nullptr;
        checked(VecGetArray(_objects.rhs, &values));
        std::copy(load.begin(), load.end(), values);
        values[0] = 0.0;
        checked(VecRestoreArray(_objects.rhs, &values));
        checked(KSPSolve(_objects.solver, _objects.rhs, _objects.solution));
        KSPConvergedReason reason;
        checked(KSPGetConvergedReason(_objects.solver, &reason));
        if (reason <= 0)
            throw std::runtime_error("Warping auxiliary solve did not converge");
        const PetscScalar* solution = nullptr;
        checked(VecGetArrayRead(_objects.solution, &solution));
        std::vector<double> result(solution, solution + static_cast<PetscInt>(_section.nodes().size()));
        checked(VecRestoreArrayRead(_objects.solution, &solution));
        double mean = 0.0;
        for (const auto& point : _section.points())
            for (std::size_t i = 0; i < 8; ++i)
                mean += point.weight * point.shape[i] * result[point.nodes[i]];
        mean /= _section.area();
        for (double& entry : result)
            entry -= mean;
        return result;
    }

    CrossSection _section;
    WarpingObjects _objects;
};

SectionWarpingSolver::SectionWarpingSolver(CrossSection section) : _impl(std::make_unique<Impl>(std::move(section))) {
}

SectionWarpingSolver::~SectionWarpingSolver() = default;

SectionWarping SectionWarpingSolver::solve(const std::vector<double>& transverse) const {
    return _impl->solve(transverse);
}

SectionAxialCorrector SectionWarpingSolver::solve_axial_corrector(const std::vector<double>& source) const {
    return _impl->solve_axial_corrector(source);
}

SectionTorsionMode build_section_torsion_mode(const CrossSection& section) {
    SectionTorsionMode result;
    result.origin = section.centroid();
    const auto n = section.nodes().size();
    result.transverse.resize(2 * n);
    for (std::size_t i = 0; i < n; ++i) {
        result.transverse[i] = -(section.nodes()[i].y - result.origin.y);
        result.transverse[n + i] = section.nodes()[i].x - result.origin.x;
    }
    result.warping = SectionWarpingSolver(section).solve(result.transverse);
    return result;
}

SectionKinematics section_torsion_kinematics(const CrossSection& section,
    const SectionTorsionMode& mode,
    std::size_t index,
    double angle,
    double twist,
    double twist_derivative) {
    const auto n = section.nodes().size();
    if (mode.transverse.size() != 2 * n || mode.warping.axial.size() != n || !std::isfinite(angle)
        || !std::isfinite(twist) || !std::isfinite(twist_derivative))
        throw std::invalid_argument("Invalid torsion mode or amplitudes");
    const auto& point = section.points().at(index);
    SectionKinematics result;
    for (std::size_t i = 0; i < 8; ++i) {
        const double ux = mode.transverse[point.nodes[i]], uy = mode.transverse[n + point.nodes[i]];
        const double axial = mode.warping.axial[point.nodes[i]];
        result.displacement[0] += point.shape[i] * ux * angle;
        result.displacement[1] += point.shape[i] * uy * angle;
        result.displacement[2] += point.shape[i] * axial * twist;
        result.gradient[0] += point.gradient[i][0] * ux * angle;
        result.gradient[1] += point.gradient[i][1] * ux * angle;
        result.gradient[2] += point.shape[i] * ux * twist;
        result.gradient[3] += point.gradient[i][0] * uy * angle;
        result.gradient[4] += point.gradient[i][1] * uy * angle;
        result.gradient[5] += point.shape[i] * uy * twist;
        result.gradient[6] += point.gradient[i][0] * axial * twist;
        result.gradient[7] += point.gradient[i][1] * axial * twist;
        result.gradient[8] += point.shape[i] * axial * twist_derivative;
    }
    return result;
}
} // namespace fuelsim
