#ifndef FUELSIM_CHECKPOINT_IO_HPP
#define FUELSIM_CHECKPOINT_IO_HPP
#include "fuelsim/transient_problem.hpp"
#include <string>
namespace fuelsim {
class TransientCheckpointIo final {
  public:
    static void write(const std::string& path, const TransientProblem& problem, double next_time_step);
    static double restore(const std::string& path, TransientProblem& problem);
};
} // namespace fuelsim
#endif
