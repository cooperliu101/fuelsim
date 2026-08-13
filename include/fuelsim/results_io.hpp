#pragma once
#include "fuelsim/mesh.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
namespace fuelsim {
UnstructuredQuad4Mesh read_exodus_quad4(const std::string& path);
void write_exodus_quad4(const std::string& path, const UnstructuredQuad4Mesh& mesh);
UnstructuredHex8Mesh read_exodus_hex8(const std::string& path);
void write_exodus_hex8(const std::string& path, const UnstructuredHex8Mesh& mesh);
void write_transient_checkpoint(const std::string& path, const TransientProblem& problem, double next_time_step);
double restore_transient_checkpoint(const std::string& path, TransientProblem& problem);
std::string next_results_segment_path(const std::string& configured_path);
class EngineeringHistoryWriter final {
  public:
    EngineeringHistoryWriter(std::string path, const TransientProblem& problem);
    void append(const TransientProblem& problem, double time_step, double next_time_step, int nonlinear_iterations);

  private:
    std::string _path;
    std::uint64_t _problem_signature;
    std::ofstream _stream;
};
void write_steady_results(const std::string& path, const UnstructuredQuad4Mesh& mesh, const SteadyProblem& problem,
    const std::vector<double>& state);
void write_steady_results(const std::string& path, const UnstructuredHex8Mesh& mesh, const SteadyProblem& problem,
    const std::vector<double>& state);
class ExodusTransientResultsWriter final {
  public:
    ExodusTransientResultsWriter(std::string path, UnstructuredQuad4Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredHex8Mesh mesh, const TransientProblem& problem);
    void append(const TransientProblem& problem);

  private:
    std::string _path;
    std::unique_ptr<UnstructuredQuad4Mesh> _rz_mesh;
    std::unique_ptr<UnstructuredHex8Mesh> _hex_mesh;
    std::uint64_t _problem_signature;
    std::size_t _step_count;
};
} // namespace fuelsim
