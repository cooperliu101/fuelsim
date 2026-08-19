#pragma once

#include <cstdint>

namespace fuelsim {
class TransientProblem;

std::uint64_t transient_problem_signature(const TransientProblem& problem);

} // namespace fuelsim
