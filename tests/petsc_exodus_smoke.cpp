#include "petsc_exodus_fixture.hpp"

#include <petscdmplex.h>

#include <cstdio>
#include <fstream>
#include <iostream>

#ifndef PETSC_HAVE_EXODUSII
#error This test requires PETSc configured with ExodusII support
#endif

namespace {

int hex_value(char digit) {
    if (digit >= '0' && digit <= '9')
        return digit - '0';
    if (digit >= 'a' && digit <= 'f')
        return digit - 'a' + 10;
    return -1;
}

bool write_fixture(const char* path) {
    const char* hex = fuelsim::test_data::two_quad_exodus_hex;
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output)
        return false;

    for (std::size_t index = 0; hex[index] != '\0'; index += 2) {
        if (hex[index + 1] == '\0')
            return false;
        const int high = hex_value(hex[index]);
        const int low = hex_value(hex[index + 1]);
        if (high < 0 || low < 0)
            return false;
        output.put(static_cast<char>((high << 4) | low));
    }
    return output.good();
}

PetscErrorCode check_quad4_mesh(DM mesh) {
    PetscInt dimension = 0;
    PetscInt cell_begin = 0;
    PetscInt cell_end = 0;
    PetscInt vertex_begin = 0;
    PetscInt vertex_end = 0;
    PetscInt cell_set_count = 0;

    PetscFunctionBeginUser;
    PetscCall(DMGetDimension(mesh, &dimension));
    PetscCall(DMPlexGetHeightStratum(mesh, 0, &cell_begin, &cell_end));
    PetscCall(DMPlexGetDepthStratum(mesh, 0, &vertex_begin, &vertex_end));
    PetscCall(DMGetLabelSize(mesh, "Cell Sets", &cell_set_count));

    PetscCheck(dimension == 2, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "Expected a 2D mesh, got dimension %" PetscInt_FMT, dimension);
    PetscCheck(cell_end - cell_begin == 2, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "Expected 2 cells, got %" PetscInt_FMT, cell_end - cell_begin);
    PetscCheck(vertex_end - vertex_begin == 6, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "Expected 6 vertices, got %" PetscInt_FMT,
               vertex_end - vertex_begin);
    PetscCheck(cell_set_count == 1, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "Expected one element block, got %" PetscInt_FMT,
               cell_set_count);

    for (PetscInt cell = cell_begin; cell < cell_end; ++cell) {
        PetscInt closure_size = 0;
        PetscInt* closure = nullptr;
        PetscInt vertex_count = 0;
        PetscInt block_id = -1;

        PetscCall(DMGetLabelValue(mesh, "Cell Sets", cell, &block_id));
        PetscCheck(block_id == 7, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                   "Cell %" PetscInt_FMT " has block ID %" PetscInt_FMT
                   " instead of 7",
                   cell, block_id);
        PetscCall(DMPlexGetTransitiveClosure(mesh, cell, PETSC_TRUE,
                                             &closure_size, &closure));
        for (PetscInt index = 0; index < closure_size; ++index) {
            const PetscInt point = closure[2 * index];
            if (point >= vertex_begin && point < vertex_end)
                ++vertex_count;
        }
        PetscCall(DMPlexRestoreTransitiveClosure(mesh, cell, PETSC_TRUE,
                                                 &closure_size, &closure));
        PetscCheck(vertex_count == 4, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                   "Cell %" PetscInt_FMT " has %" PetscInt_FMT
                   " vertices instead of 4",
                   cell, vertex_count);
    }

    PetscCall(PetscPrintf(
        PETSC_COMM_WORLD,
        "Exodus Quad4 smoke: dimension=%" PetscInt_FMT " cells=%" PetscInt_FMT
        " vertices=%" PetscInt_FMT " block=7\n",
        dimension, cell_end - cell_begin, vertex_end - vertex_begin));
    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode run_smoke(const char* path) {
    DM loaded_mesh = nullptr;

    PetscFunctionBeginUser;
    PetscCheck(write_fixture(path), PETSC_COMM_WORLD, PETSC_ERR_FILE_WRITE,
               "Could not write the independent Exodus fixture");
    PetscCall(DMPlexCreateExodusFromFile(PETSC_COMM_WORLD, path, PETSC_TRUE,
                                         &loaded_mesh));
    PetscCall(check_quad4_mesh(loaded_mesh));
    PetscCall(DMDestroy(&loaded_mesh));
    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_petsc_exodus_smoke <output.exo>\n";
        return 2;
    }

    PetscErrorCode error = PetscInitialize(&argc, &argv, nullptr, nullptr);
    if (error != PETSC_SUCCESS)
        return static_cast<int>(error);

    error = run_smoke(argv[1]);
    const PetscErrorCode finalize_error = PetscFinalize();
    const int remove_error = std::remove(argv[1]);
    if (error != PETSC_SUCCESS)
        return static_cast<int>(error);
    if (finalize_error != PETSC_SUCCESS)
        return static_cast<int>(finalize_error);
    if (remove_error != 0) {
        std::cerr << "Failed to remove smoke-test Exodus file\n";
        return 3;
    }
    return 0;
}
