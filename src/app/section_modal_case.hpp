#pragma once
#include "io/case_input.hpp"
#include "solver/petsc_solver.hpp"

namespace fuelsim {
void run_section_modal_case(const FuelSimCaseDefinition& definition, const PetscSession& session);
}
