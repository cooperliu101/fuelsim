#include "modal_condensation.hpp"
#include <stdexcept>
#include <string>

namespace fuelsim {
namespace {
void check(PetscErrorCode code) {
    if (code != PETSC_SUCCESS)
        throw std::runtime_error("Modal condensation PETSc failure: " + std::to_string(code));
}

struct Objects final {
    IS local_ids = nullptr, retained_ids = nullptr;
    Mat local = nullptr, coupling = nullptr, schur = nullptr, responses = nullptr, response_rhs = nullptr;
    KSP local_solver = nullptr, retained_solver = nullptr;
    Vec local_rhs = nullptr, local_solution = nullptr;
    Vec retained_rhs = nullptr, retained_solution = nullptr, retained_work = nullptr;

    ~Objects() {
        (void)MatDestroy(&responses);
        (void)MatDestroy(&response_rhs);
        (void)VecDestroy(&retained_work);
        (void)VecDestroy(&retained_solution);
        (void)VecDestroy(&retained_rhs);
        (void)VecDestroy(&local_solution);
        (void)VecDestroy(&local_rhs);
        (void)KSPDestroy(&retained_solver);
        (void)KSPDestroy(&local_solver);
        (void)MatDestroy(&schur);
        (void)MatDestroy(&coupling);
        (void)MatDestroy(&local);
        (void)ISDestroy(&retained_ids);
        (void)ISDestroy(&local_ids);
    }
};

std::vector<double> gather(Vec vector) {
    Vec all = nullptr;
    VecScatter scatter = nullptr;
    check(VecScatterCreateToAll(vector, &scatter, &all));
    check(VecScatterBegin(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD));
    check(VecScatterEnd(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD));
    PetscInt size;
    check(VecGetSize(all, &size));
    const PetscScalar* values;
    check(VecGetArrayRead(all, &values));
    std::vector<double> result(values, values + size);
    check(VecRestoreArrayRead(all, &values));
    check(VecScatterDestroy(&scatter));
    check(VecDestroy(&all));
    return result;
}

void direct_solver(Mat matrix, KSP* solver) {
    check(KSPCreate(PETSC_COMM_WORLD, solver));
    check(KSPSetOperators(*solver, matrix, matrix));
    check(KSPSetType(*solver, KSPPREONLY));
    check(KSPSetErrorIfNotConverged(*solver, PETSC_TRUE));
    PC pc;
    check(KSPGetPC(*solver, &pc));
    check(PCSetType(pc, PCLU));
    check(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
}

void solve_direct(KSP solver, Vec rhs, Vec solution) {
    check(KSPSolve(solver, rhs, solution));
    KSPConvergedReason reason;
    check(KSPGetConvergedReason(solver, &reason));
    if (reason <= 0)
        throw std::runtime_error("Modal local or condensed factorization failed");
}
} // namespace

class ModalStaticCondensation::Impl final {
  public:
    Impl(Mat matrix, const std::vector<bool>& eliminated) {
        PetscInt size, begin, end;
        check(MatGetSize(matrix, &size, nullptr));
        check(MatGetOwnershipRange(matrix, &begin, &end));
        if (eliminated.size() != static_cast<std::size_t>(size))
            throw std::invalid_argument("Condensation partition has the wrong size");
        std::vector<PetscInt> local, retained;
        std::size_t local_count = 0;
        for (PetscInt i = 0; i < size; ++i) {
            if (eliminated[static_cast<std::size_t>(i)])
                ++local_count;
            if (i >= begin && i < end)
                (eliminated[static_cast<std::size_t>(i)] ? local : retained).push_back(i);
        }
        if (local_count == 0 || local_count == eliminated.size())
            throw std::invalid_argument("Condensation requires both local and retained unknowns");
        check(ISCreateGeneral(PETSC_COMM_WORLD,
            static_cast<PetscInt>(local.size()),
            local.data(),
            PETSC_COPY_VALUES,
            &_objects.local_ids));
        check(ISCreateGeneral(PETSC_COMM_WORLD,
            static_cast<PetscInt>(retained.size()),
            retained.data(),
            PETSC_COPY_VALUES,
            &_objects.retained_ids));
        check(MatCreateSubMatrix(matrix, _objects.local_ids, _objects.local_ids, MAT_INITIAL_MATRIX, &_objects.local));
        check(MatCreateSubMatrix(matrix,
            _objects.local_ids,
            _objects.retained_ids,
            MAT_INITIAL_MATRIX,
            &_objects.coupling));
        check(MatCreateSubMatrix(matrix,
            _objects.retained_ids,
            _objects.retained_ids,
            MAT_INITIAL_MATRIX,
            &_objects.schur));
        check(MatSetOption(_objects.schur, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
        check(MatCreateVecs(_objects.local, &_objects.local_solution, &_objects.local_rhs));
        check(MatCreateVecs(_objects.schur, &_objects.retained_solution, &_objects.retained_rhs));
        check(VecDuplicate(_objects.retained_rhs, &_objects.retained_work));
        direct_solver(_objects.local, &_objects.local_solver);

        // Only columns touching an end region need a local response solve.
        // Keep every nonzero coupling, including small entries; no numerical
        // threshold is used to reduce the interface or alter its stiffness.
        check(VecSet(_objects.retained_work, 0.0));
        check(MatGetOwnershipRange(_objects.coupling, &begin, &end));
        for (PetscInt i = begin; i < end; ++i) {
            PetscInt count;
            const PetscInt* columns;
            const PetscScalar* values;
            check(MatGetRow(_objects.coupling, i, &count, &columns, &values));
            for (PetscInt j = 0; j < count; ++j)
                if (values[j] != 0.0)
                    check(VecSetValue(_objects.retained_work, columns[j], 1.0, ADD_VALUES));
            check(MatRestoreRow(_objects.coupling, i, &count, &columns, &values));
        }
        check(VecAssemblyBegin(_objects.retained_work));
        check(VecAssemblyEnd(_objects.retained_work));
        const auto flags = gather(_objects.retained_work);
        for (std::size_t i = 0; i < flags.size(); ++i)
            if (flags[i] != 0.0)
                _ports.push_back(static_cast<PetscInt>(i));
        if (_ports.empty())
            throw std::invalid_argument("Local modal region has no retained interface");
        // All interface responses share one factorization. Supply them as a
        // distributed dense RHS block so MUMPS can traverse its factors once
        // for multiple columns, rather than issuing one scalar solve per port.
        PetscInt local_rows, rows;
        check(MatGetLocalSize(_objects.local, &local_rows, nullptr));
        check(MatGetSize(_objects.local, &rows, nullptr));
        check(MatCreateDense(PETSC_COMM_WORLD,
            local_rows,
            PETSC_DECIDE,
            rows,
            static_cast<PetscInt>(_ports.size()),
            nullptr,
            &_objects.response_rhs));
        check(MatZeroEntries(_objects.response_rhs));
        check(MatDuplicate(_objects.response_rhs, MAT_DO_NOT_COPY_VALUES, &_objects.responses));
        PetscScalar* columns;
        PetscInt leading;
        check(MatDenseGetLDA(_objects.response_rhs, &leading));
        check(MatDenseGetArray(_objects.response_rhs, &columns));
        std::vector<PetscInt> port_index(flags.size(), -1);
        for (std::size_t j = 0; j < _ports.size(); ++j)
            port_index[static_cast<std::size_t>(_ports[j])] = static_cast<PetscInt>(j);
        check(MatGetOwnershipRange(_objects.coupling, &begin, &end));
        for (PetscInt i = begin; i < end; ++i) {
            PetscInt count;
            const PetscInt* ids;
            const PetscScalar* values;
            check(MatGetRow(_objects.coupling, i, &count, &ids, &values));
            for (PetscInt j = 0; j < count; ++j)
                if (values[j] != 0.0)
                    columns[port_index[static_cast<std::size_t>(ids[j])] * leading + i - begin] = values[j];
            check(MatRestoreRow(_objects.coupling, i, &count, &ids, &values));
        }
        check(MatDenseRestoreArray(_objects.response_rhs, &columns));
        check(MatAssemblyBegin(_objects.response_rhs, MAT_FINAL_ASSEMBLY));
        check(MatAssemblyEnd(_objects.response_rhs, MAT_FINAL_ASSEMBLY));
        check(KSPMatSolve(_objects.local_solver, _objects.response_rhs, _objects.responses));
        check(MatDestroy(&_objects.response_rhs));
        check(MatGetOwnershipRange(_objects.schur, &begin, &end));
        for (std::size_t j = 0; j < _ports.size(); ++j) {
            const auto port = _ports[j];
            Vec response = nullptr;
            check(MatDenseGetColumnVecRead(_objects.responses, static_cast<PetscInt>(j), &response));
            check(MatMultTranspose(_objects.coupling, response, _objects.retained_work));
            check(MatDenseRestoreColumnVecRead(_objects.responses, static_cast<PetscInt>(j), &response));
            const PetscScalar* values;
            check(VecGetArrayRead(_objects.retained_work, &values));
            for (PetscInt i = begin; i < end; ++i)
                if (values[i - begin] != 0.0)
                    check(MatSetValue(_objects.schur, i, port, -values[i - begin], ADD_VALUES));
            check(VecRestoreArrayRead(_objects.retained_work, &values));
        }
        check(MatAssemblyBegin(_objects.schur, MAT_FINAL_ASSEMBLY));
        check(MatAssemblyEnd(_objects.schur, MAT_FINAL_ASSEMBLY));
        direct_solver(_objects.schur, &_objects.retained_solver);
    }

    void solve(Vec rhs, Vec solution) const {
        Vec view = nullptr;
        check(VecGetSubVector(rhs, _objects.local_ids, &view));
        check(VecCopy(view, _objects.local_rhs));
        check(VecRestoreSubVector(rhs, _objects.local_ids, &view));
        check(VecGetSubVector(rhs, _objects.retained_ids, &view));
        check(VecCopy(view, _objects.retained_rhs));
        check(VecRestoreSubVector(rhs, _objects.retained_ids, &view));
        solve_direct(_objects.local_solver, _objects.local_rhs, _objects.local_solution);
        check(MatMultTranspose(_objects.coupling, _objects.local_solution, _objects.retained_work));
        check(VecAXPY(_objects.retained_rhs, -1.0, _objects.retained_work));
        solve_direct(_objects.retained_solver, _objects.retained_rhs, _objects.retained_solution);
        const auto retained = gather(_objects.retained_solution);
        for (std::size_t j = 0; j < _ports.size(); ++j) {
            Vec response = nullptr;
            check(MatDenseGetColumnVecRead(_objects.responses, static_cast<PetscInt>(j), &response));
            check(VecAXPY(_objects.local_solution, -retained[static_cast<std::size_t>(_ports[j])], response));
            check(MatDenseRestoreColumnVecRead(_objects.responses, static_cast<PetscInt>(j), &response));
        }
        check(VecGetSubVector(solution, _objects.local_ids, &view));
        check(VecCopy(_objects.local_solution, view));
        check(VecRestoreSubVector(solution, _objects.local_ids, &view));
        check(VecGetSubVector(solution, _objects.retained_ids, &view));
        check(VecCopy(_objects.retained_solution, view));
        check(VecRestoreSubVector(solution, _objects.retained_ids, &view));
    }

  private:
    Objects _objects;
    std::vector<PetscInt> _ports;
};

ModalStaticCondensation::ModalStaticCondensation(Mat matrix, const std::vector<bool>& eliminated)
    : _impl(std::make_unique<Impl>(matrix, eliminated)) {
}

ModalStaticCondensation::~ModalStaticCondensation() = default;

void ModalStaticCondensation::solve(Vec rhs, Vec solution) const {
    _impl->solve(rhs, solution);
}
} // namespace fuelsim
