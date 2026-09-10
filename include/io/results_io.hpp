#pragma once
#include "core/mesh.hpp"
#include "core/steady_problem.hpp"
#include "core/transient_problem.hpp"
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>

namespace fuelsim {
UnstructuredQuad4Mesh read_exodus_quad4(const std::string& path);
void write_exodus_quad4(const std::string& path, const UnstructuredQuad4Mesh& mesh);
UnstructuredQuad8Mesh read_exodus_quad8(const std::string& path);
void write_exodus_quad8(const std::string& path, const UnstructuredQuad8Mesh& mesh);
bool exodus_uses_quad8(const std::string& path);
UnstructuredHex8Mesh read_exodus_hex8(const std::string& path);
void write_exodus_hex8(const std::string& path, const UnstructuredHex8Mesh& mesh);
UnstructuredHex20Mesh read_exodus_hex20(const std::string& path);
void write_exodus_hex20(const std::string& path, const UnstructuredHex20Mesh& mesh);
bool exodus_uses_hex20(const std::string& path);
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

void write_steady_results(const std::string& path,
    const UnstructuredQuad4Mesh& mesh,
    const SteadyProblem& problem,
    const std::vector<double>& state);
void write_steady_results(const std::string& path,
    const UnstructuredQuad8Mesh& mesh,
    const SteadyProblem& problem,
    const std::vector<double>& state);
void write_steady_results(const std::string& path,
    const UnstructuredHex8Mesh& mesh,
    const SteadyProblem& problem,
    const std::vector<double>& state);
void write_steady_results(const std::string& path,
    const UnstructuredHex20Mesh& mesh,
    const SteadyProblem& problem,
    const std::vector<double>& state);

class ExodusTransientResultsWriter final {
  public:
    ExodusTransientResultsWriter(std::string path, UnstructuredQuad4Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredQuad8Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredHex8Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredHex20Mesh mesh, const TransientProblem& problem);
    void append(const TransientProblem& problem);

  private:
    std::string _path;
    std::unique_ptr<UnstructuredQuad4Mesh> _rz_mesh;
    std::unique_ptr<UnstructuredQuad8Mesh> _quad8_mesh;
    std::unique_ptr<UnstructuredHex8Mesh> _hex_mesh;
    std::unique_ptr<UnstructuredHex20Mesh> _hex20_mesh;
    std::uint64_t _problem_signature;
    std::size_t _step_count;
};
} // namespace fuelsim
