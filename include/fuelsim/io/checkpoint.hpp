#pragma once

#include "fuelsim/core/mesh.hpp"
#include "fuelsim/core/transient_problem.hpp"

#include <string>

namespace fuelsim {

void write_transient_checkpoint(const std::string& path, const TransientProblem& problem, double next_time_step);
double restore_transient_checkpoint(const std::string& path, TransientProblem& problem);

} // namespace fuelsim
