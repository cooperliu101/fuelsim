#ifndef FUELSIM_TRANSIENT_PROBLEM_SIGNATURE_HPP
#define FUELSIM_TRANSIENT_PROBLEM_SIGNATURE_HPP

#include <cstdint>

namespace fuelsim {

class TransientProblem;

std::uint64_t transient_problem_signature(const TransientProblem& problem);

} // namespace fuelsim

#endif
