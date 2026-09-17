#pragma once
#include "core/mesh.hpp"
#include "core/plane_mesh.hpp"
#include "core/steady_problem.hpp"
#include "core/transient_problem.hpp"
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>

namespace fuelsim {
UnstructuredBar2Mesh read_exodus_bar2(const std::string& path);
void write_exodus_bar2(const std::string& path, const UnstructuredBar2Mesh& mesh);
UnstructuredQuad4Mesh read_exodus_quad4(const std::string& path);
void write_exodus_quad4(const std::string& path, const UnstructuredQuad4Mesh& mesh);
UnstructuredPlaneQuad8Mesh read_exodus_plane_quad8(const std::string& path);
void write_exodus_plane_quad8(const std::string& path, const UnstructuredPlaneQuad8Mesh& mesh);
void write_steady_plane_results(const std::string& path,
    const UnstructuredPlaneQuad8Mesh& mesh,
    const SteadyProblem& problem,
    const std::vector<double>& state);
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
    const UnstructuredBar2Mesh& mesh,
    const SteadyProblem& problem,
    const std::vector<double>& state);
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
    ExodusTransientResultsWriter(std::string path, UnstructuredPlaneQuad8Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredBar2Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredQuad4Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredQuad8Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredHex8Mesh mesh, const TransientProblem& problem);
    ExodusTransientResultsWriter(std::string path, UnstructuredHex20Mesh mesh, const TransientProblem& problem);
    void append(const TransientProblem& problem);

  private:
    std::string _path;
    std::unique_ptr<UnstructuredPlaneQuad8Mesh> _plane_mesh;
    std::unique_ptr<UnstructuredBar2Mesh> _bar2_mesh;
    std::unique_ptr<UnstructuredQuad4Mesh> _rz_mesh;
    std::unique_ptr<UnstructuredQuad8Mesh> _quad8_mesh;
    std::unique_ptr<UnstructuredHex8Mesh> _hex_mesh;
    std::unique_ptr<UnstructuredHex20Mesh> _hex20_mesh;
    std::uint64_t _problem_signature;
    std::size_t _step_count;
};
} // namespace fuelsim
