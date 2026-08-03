#ifndef FUELSIM_RESULTS_IO_HPP
#define FUELSIM_RESULTS_IO_HPP

#include "fuelsim/mesh.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>

namespace fuelsim {

std::string next_results_segment_path(const std::string& configured_path);

class EngineeringHistoryWriter final {
  public:
    EngineeringHistoryWriter(std::string path,
                             const TransientProblem& problem);

    void append(const TransientProblem& problem, double time_step,
                double next_time_step, int nonlinear_iterations);

  private:
    std::string _path;
    std::uint64_t _problem_signature;
    std::ofstream _stream;
};

class ExodusResultsIo final {
  public:
    static void write_steady(const std::string& path,
                             const UnstructuredQuad4Mesh& mesh,
                             const SteadyProblem& problem,
                             const std::vector<double>& state);
};

class ExodusTransientResultsWriter final {
  public:
    ExodusTransientResultsWriter(std::string path, UnstructuredQuad4Mesh mesh,
                                 const TransientProblem& problem);

    void append(const TransientProblem& problem);
    std::size_t step_count() const noexcept;

  private:
    std::string _path;
    UnstructuredQuad4Mesh _mesh;
    std::uint64_t _problem_signature;
    std::size_t _step_count;
};

} // namespace fuelsim

#endif
