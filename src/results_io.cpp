#include "fuelsim/results_io.hpp"

#include "fuelsim/exodus_mesh_io.hpp"

#include <exodusII.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
namespace {

void check_exodus(int status, const std::string& operation) {
    if (status < 0)
        throw std::runtime_error(operation + ": " + ex_strerror(status));
}

int checked_int(std::size_t value, const std::string& quantity) {
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error(quantity + " exceeds the Exodus int range");
    return static_cast<int>(value);
}

class ExodusFile final {
  public:
    explicit ExodusFile(int id) : _id(id) {}
    ExodusFile(const ExodusFile&) = delete;
    ExodusFile& operator=(const ExodusFile&) = delete;

    ~ExodusFile() {
        if (_id >= 0)
            ex_close(_id);
    }

    int id() const noexcept {
        return _id;
    }

    void close() {
        const int id = _id;
        _id = -1;
        check_exodus(ex_close(id), "Could not close Exodus results file");
    }

  private:
    int _id;
};

ExodusFile open_results(const std::string& path) {
    int cpu_word_size = static_cast<int>(sizeof(double));
    int io_word_size = 0;
    float version = 0.0F;
    const int exoid = ex_open(path.c_str(), EX_WRITE, &cpu_word_size,
                              &io_word_size, &version);
    if (exoid < 0)
        throw std::runtime_error("Could not open Exodus results file '" + path +
                                 "': " + ex_strerror(exoid));
    ex_set_int64_status(exoid, EX_ALL_INT64_API);
    return ExodusFile(exoid);
}

std::vector<char*> variable_name_pointers(std::vector<std::string>& names) {
    std::vector<char*> result;
    result.reserve(names.size());
    for (std::string& name : names)
        result.push_back(name.data());
    return result;
}

std::vector<std::string>
nodal_variable_names(const std::vector<ContactDefinition>& contacts) {
    std::vector<std::string> result = {"temperature", "displacement_r",
                                       "displacement_z"};
    for (const ContactDefinition& contact : contacts) {
        result.push_back("contact_gap_" + contact.name);
        result.push_back("contact_pressure_" + contact.name);
        result.push_back("contact_tangential_traction_" + contact.name);
        result.push_back("contact_elastic_tangential_slip_" + contact.name);
        result.push_back("contact_sliding_" + contact.name);
    }
    return result;
}

std::vector<std::string>
global_variable_names(const std::vector<ContactDefinition>& contacts) {
    std::vector<std::string> result = {"load_factor"};
    for (const ContactDefinition& contact : contacts) {
        result.push_back("contact_heat_rate_" + contact.name);
        result.push_back("contact_force_" + contact.name);
        result.push_back("contact_tangential_force_" + contact.name);
    }
    return result;
}

std::vector<std::string> stress_variable_names() {
    const std::array<const char*, 4> components = {"rr", "zz", "hoop", "rz"};
    std::vector<std::string> result;
    result.reserve(16);
    for (std::size_t q = 0; q < 4; ++q) {
        for (const char* component : components)
            result.push_back("stress_" + std::string(component) + "_q" +
                             std::to_string(q));
    }
    return result;
}

std::vector<std::string> transient_element_variable_names() {
    const std::array<const char*, 4> components = {"rr", "zz", "hoop", "rz"};
    std::vector<std::string> result = stress_variable_names();
    result.reserve(56);
    for (std::size_t q = 0; q < 4; ++q) {
        for (const char* component : components)
            result.push_back("plastic_" + std::string(component) + "_q" +
                             std::to_string(q));
        for (const char* component : components)
            result.push_back("creep_" + std::string(component) + "_q" +
                             std::to_string(q));
        result.push_back("equiv_plastic_q" + std::to_string(q));
        result.push_back("equiv_creep_q" + std::to_string(q));
    }
    return result;
}

void define_variables(const std::string& path,
                      const UnstructuredQuad4Mesh& mesh,
                      std::vector<std::string> nodal_names,
                      std::vector<std::string> element_names,
                      std::vector<std::string> global_names) {
    ExodusMeshIo::write_quad4(path, mesh);
    ExodusFile file = open_results(path);
    check_exodus(ex_set_max_name_length(file.id(), 64),
                 "Could not set Exodus result-name length");

    std::vector<char*> nodal = variable_name_pointers(nodal_names);
    check_exodus(ex_put_variable_param(file.id(), EX_NODAL,
                                       static_cast<int>(nodal.size())),
                 "Could not define Exodus nodal variables");
    check_exodus(ex_put_variable_names(file.id(), EX_NODAL,
                                       static_cast<int>(nodal.size()),
                                       nodal.data()),
                 "Could not name Exodus nodal variables");

    std::vector<char*> global = variable_name_pointers(global_names);
    check_exodus(ex_put_variable_param(file.id(), EX_GLOBAL,
                                       static_cast<int>(global.size())),
                 "Could not define Exodus global variables");
    check_exodus(ex_put_variable_names(file.id(), EX_GLOBAL,
                                       static_cast<int>(global.size()),
                                       global.data()),
                 "Could not name Exodus global variables");

    std::vector<char*> element = variable_name_pointers(element_names);
    check_exodus(ex_put_variable_param(file.id(), EX_ELEM_BLOCK,
                                       static_cast<int>(element.size())),
                 "Could not define Exodus element variables");
    check_exodus(ex_put_variable_names(file.id(), EX_ELEM_BLOCK,
                                       static_cast<int>(element.size()),
                                       element.data()),
                 "Could not name Exodus element variables");
    std::vector<int> truth(mesh.element_blocks().size() * element.size(), 1);
    check_exodus(
        ex_put_truth_table(file.id(), EX_ELEM_BLOCK,
                           static_cast<int>(mesh.element_blocks().size()),
                           static_cast<int>(element.size()), truth.data()),
        "Could not define Exodus element-variable truth table");
    file.close();
}

std::vector<std::size_t> block_elements(const UnstructuredQuad4Mesh& mesh,
                                        std::int64_t block_id) {
    std::vector<std::size_t> result;
    for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
        if (mesh.element_block_ids()[element] == block_id)
            result.push_back(element);
    }
    return result;
}

void write_step(const std::string& path, const UnstructuredQuad4Mesh& mesh,
                std::size_t step, double time,
                const std::vector<std::vector<double>>& nodal_values,
                const std::vector<std::vector<double>>& element_values,
                const std::vector<double>& global_values) {
    if (step == 0 || nodal_values.empty() || element_values.empty() ||
        global_values.empty())
        throw std::invalid_argument("Exodus result step is incomplete");
    ExodusFile file = open_results(path);
    const int exodus_step = checked_int(step, "Exodus result step");
    check_exodus(ex_put_time(file.id(), exodus_step, &time),
                 "Could not write Exodus result time");
    check_exodus(ex_put_var(file.id(), exodus_step, EX_GLOBAL, 1, 0,
                            static_cast<std::int64_t>(global_values.size()),
                            global_values.data()),
                 "Could not write Exodus global results");
    for (std::size_t variable = 0; variable < nodal_values.size(); ++variable) {
        if (nodal_values[variable].size() != mesh.nodes().size())
            throw std::invalid_argument(
                "Exodus nodal result size does not match mesh");
        check_exodus(ex_put_var(file.id(), exodus_step, EX_NODAL,
                                static_cast<int>(variable + 1), 1,
                                static_cast<std::int64_t>(mesh.nodes().size()),
                                nodal_values[variable].data()),
                     "Could not write Exodus nodal results");
    }
    for (std::size_t variable = 0; variable < element_values.size();
         ++variable) {
        if (element_values[variable].size() != mesh.elements().size())
            throw std::invalid_argument(
                "Exodus element result size does not match mesh");
        for (const ElementBlockInfo& block : mesh.element_blocks()) {
            const std::vector<std::size_t> elements =
                block_elements(mesh, block.id);
            std::vector<double> values;
            values.reserve(elements.size());
            for (const std::size_t element : elements)
                values.push_back(element_values[variable][element]);
            check_exodus(ex_put_var(file.id(), exodus_step, EX_ELEM_BLOCK,
                                    static_cast<int>(variable + 1), block.id,
                                    static_cast<std::int64_t>(values.size()),
                                    values.data()),
                         "Could not write Exodus element results");
        }
    }
    check_exodus(ex_update(file.id()), "Could not flush Exodus results");
    file.close();
}

void fill_steady_nodal(const UnstructuredQuad4Mesh& mesh,
                       const SteadyProblem& problem,
                       const std::vector<double>& state,
                       std::vector<std::vector<double>>& values) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    values.assign(nodal_variable_names(problem.definition().contacts).size(),
                  std::vector<double>(mesh.nodes().size(), missing));
    std::vector<bool> present(mesh.nodes().size(), false);
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const RegionMesh& region_mesh = problem.region_mesh(region);
        const std::size_t offset = problem.region_node_offset(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size();
             ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            if (present.at(source))
                throw std::invalid_argument(
                    "Exodus result mapping contains a shared source node");
            present[source] = true;
            const std::size_t global = offset + local;
            values[0][source] = state.at(problem.dof_map().temperature(global));
            values[1][source] =
                state.at(problem.dof_map().radial_displacement(global));
            values[2][source] =
                state.at(problem.dof_map().axial_displacement(global));
        }
    }
    for (std::size_t contact = 0; contact < problem.contact_count();
         ++contact) {
        const std::vector<std::size_t> nodes =
            problem.contact_secondary_source_nodes(contact);
        const std::vector<ContactNodeSummary> summary =
            problem.summarize_contact_nodes(contact, state);
        if (nodes.size() != summary.size())
            throw std::logic_error("Contact result mapping size mismatch");
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            const std::size_t base = 3 + 5 * contact;
            values[base].at(nodes[node]) =
                summary[node].projected ? summary[node].gap : missing;
            values[base + 1].at(nodes[node]) =
                summary[node].projected ? summary[node].pressure : missing;
            values[base + 2].at(nodes[node]) =
                summary[node].projected ? summary[node].tangential_traction
                                        : missing;
            values[base + 3].at(nodes[node]) =
                summary[node].projected
                    ? summary[node].elastic_tangential_slip
                    : missing;
            values[base + 4].at(nodes[node]) =
                summary[node].projected
                    ? (summary[node].sliding ? 1.0 : 0.0)
                    : missing;
        }
    }
}

void fill_transient_nodal(const UnstructuredQuad4Mesh& mesh,
                          const TransientProblem& problem,
                          std::vector<std::vector<double>>& values) {
    const std::vector<ContactDefinition>& contact_definitions =
        problem.definition().spatial.contacts;
    const std::size_t contacts = contact_definitions.size();
    const double missing = std::numeric_limits<double>::quiet_NaN();
    values.assign(nodal_variable_names(contact_definitions).size(),
                  std::vector<double>(mesh.nodes().size(), missing));
    std::vector<bool> present(mesh.nodes().size(), false);
    const std::vector<double>& state = problem.committed_solution();
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const RegionMesh& region_mesh = problem.region_mesh(region);
        const std::size_t offset = problem.region_node_offset(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size();
             ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            if (present.at(source))
                throw std::invalid_argument(
                    "Exodus result mapping contains a shared source node");
            present[source] = true;
            const std::size_t global = offset + local;
            values[0][source] = state.at(problem.dof_map().temperature(global));
            values[1][source] =
                state.at(problem.dof_map().radial_displacement(global));
            values[2][source] =
                state.at(problem.dof_map().axial_displacement(global));
        }
    }
    for (std::size_t contact = 0; contact < contacts; ++contact) {
        const std::vector<std::size_t> nodes =
            problem.contact_secondary_source_nodes(contact);
        const std::vector<ContactNodeSummary> summary =
            problem.summarize_contact_nodes(contact, state);
        if (nodes.size() != summary.size())
            throw std::logic_error("Contact result mapping size mismatch");
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            const std::size_t base = 3 + 5 * contact;
            values[base].at(nodes[node]) =
                summary[node].projected ? summary[node].gap : missing;
            values[base + 1].at(nodes[node]) =
                summary[node].projected ? summary[node].pressure : missing;
            values[base + 2].at(nodes[node]) =
                summary[node].projected ? summary[node].tangential_traction
                                        : missing;
            values[base + 3].at(nodes[node]) =
                summary[node].projected
                    ? summary[node].elastic_tangential_slip
                    : missing;
            values[base + 4].at(nodes[node]) =
                summary[node].projected
                    ? (summary[node].sliding ? 1.0 : 0.0)
                    : missing;
        }
    }
}

std::vector<double> steady_globals(const SteadyProblem& problem,
                                   const std::vector<double>& state) {
    std::vector<double> result = {problem.load_factor()};
    for (std::size_t contact = 0; contact < problem.contact_count();
         ++contact) {
        const InterfaceSummary summary =
            problem.summarize_interface(contact, state);
        result.push_back(summary.total_heat_rate);
        result.push_back(summary.total_contact_force);
        result.push_back(summary.total_tangential_force);
    }
    return result;
}

std::vector<double> transient_globals(const TransientProblem& problem) {
    std::vector<double> result = {problem.committed_load_factor()};
    const std::vector<double>& state = problem.committed_solution();
    const std::size_t contacts = problem.definition().spatial.contacts.size();
    for (std::size_t contact = 0; contact < contacts; ++contact) {
        const InterfaceSummary summary =
            problem.summarize_interface(contact, state);
        result.push_back(summary.total_heat_rate);
        result.push_back(summary.total_contact_force);
        result.push_back(summary.total_tangential_force);
    }
    return result;
}

std::vector<std::vector<double>>
steady_elements(const UnstructuredQuad4Mesh& mesh, const SteadyProblem& problem,
                const std::vector<double>& state) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> result(
        stress_variable_names().size(),
        std::vector<double>(mesh.elements().size(), missing));
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const RegionMesh& region_mesh = problem.region_mesh(region);
        const std::size_t contribution_offset =
            problem.region_element_offset(region);
        for (std::size_t element = 0; element < region_mesh.elements().size();
             ++element) {
            const LocalValues local = problem.contribution_state(
                contribution_offset + element, state);
            const auto stresses = problem.region_kernel(region).stress_values(
                problem.region_element_geometry(region, element), local);
            const std::size_t source =
                region_mesh.source_element_ids().at(element);
            for (std::size_t q = 0; q < 4; ++q) {
                result[4 * q][source] = stresses[q].rr;
                result[4 * q + 1][source] = stresses[q].zz;
                result[4 * q + 2][source] = stresses[q].hoop;
                result[4 * q + 3][source] = stresses[q].rz;
            }
        }
    }
    return result;
}

std::vector<std::vector<double>>
transient_elements(const UnstructuredQuad4Mesh& mesh,
                   const TransientProblem& problem) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> result(
        transient_element_variable_names().size(),
        std::vector<double>(mesh.elements().size(), missing));
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const RegionMesh& region_mesh = problem.region_mesh(region);
        for (std::size_t element = 0; element < region_mesh.elements().size();
             ++element) {
            const std::size_t source =
                region_mesh.source_element_ids().at(element);
            const Quad4MaterialHistory& history =
                problem.material_history(region, element);
            const auto& stresses = problem.material_stress(region, element);
            for (std::size_t q = 0; q < 4; ++q) {
                const std::size_t stress = 4 * q;
                result[stress][source] = stresses[q].rr;
                result[stress + 1][source] = stresses[q].zz;
                result[stress + 2][source] = stresses[q].hoop;
                result[stress + 3][source] = stresses[q].rz;
                const std::size_t history_offset = 16 + 10 * q;
                for (std::size_t component = 0; component < 4; ++component) {
                    result[history_offset + component][source] =
                        history[q].plastic_strain[component];
                    result[history_offset + 4 + component][source] =
                        history[q].creep_strain[component];
                }
                result[history_offset + 8][source] =
                    history[q].equivalent_plastic_strain;
                result[history_offset + 9][source] =
                    history[q].equivalent_creep_strain;
            }
        }
    }
    return result;
}

} // namespace

std::string next_results_segment_path(const std::string& configured_path) {
    if (configured_path.empty())
        throw std::invalid_argument(
            "Results segment path requires a configured path");
    const std::filesystem::path configured(configured_path);
    const std::filesystem::path directory = configured.parent_path();
    const std::string stem = configured.stem().string();
    const std::string extension = configured.extension().string();
    for (std::size_t segment = 1;; ++segment) {
        const std::filesystem::path candidate =
            directory /
            (stem + ".part" + std::to_string(segment) + extension);
        if (!std::filesystem::exists(candidate))
            return candidate.string();
    }
}

EngineeringHistoryWriter::EngineeringHistoryWriter(
    std::string path, const TransientProblem& problem)
    : _path(std::move(path)),
      _problem_signature(problem.committed_state_signature()) {
    if (_path.empty())
        throw std::invalid_argument(
            "Engineering history path must not be empty");
    _stream.open(_path, std::ios::out | std::ios::trunc);
    if (!_stream)
        throw std::runtime_error(
            "Could not open engineering history file '" + _path + "'");
    _stream.exceptions(std::ios::badbit | std::ios::failbit);
    _stream << "time,time_step,next_time_step,load_factor,"
               "nonlinear_iterations,generated_heat_rate,stored_heat_rate,"
               "convection_heat_rate,interface_heat_imbalance,"
               "dirichlet_heat_input_rate,global_thermal_balance,"
               "relative_thermal_balance,unconstrained_thermal_residual_l2,"
               "internal_mechanical_work_increment,"
               "pressure_traction_work_increment,"
               "dirichlet_reaction_work_increment,contact_work_increment,"
               "mechanical_work_balance,relative_mechanical_work_balance,"
               "unconstrained_mechanical_residual_l2,elastic_energy_change,"
               "plastic_dissipation_increment,creep_dissipation_increment";
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const std::string prefix = ",region_" + problem.region(region).name;
        _stream << prefix << "_maximum_temperature" << prefix
                << "_maximum_equivalent_plastic_strain" << prefix
                << "_maximum_equivalent_creep_strain";
    }
    for (const ContactDefinition& contact :
         problem.definition().spatial.contacts) {
        const std::string prefix = ",contact_" + contact.name;
        _stream << prefix << "_minimum_gap" << prefix
                << "_maximum_pressure" << prefix << "_total_heat_rate"
                << prefix << "_total_force" << prefix
                << "_total_tangential_force";
    }
    _stream << '\n' << std::scientific << std::setprecision(12);
}

void EngineeringHistoryWriter::append(const TransientProblem& problem,
                                      double time_step,
                                      double next_time_step,
                                      int nonlinear_iterations) {
    if (problem.time_step_active())
        throw std::logic_error(
            "Engineering history cannot be written during an active time step");
    if (problem.committed_state_signature() != _problem_signature)
        throw std::invalid_argument(
            "Engineering history problem does not match writer model");
    _stream << problem.committed_time() << ',' << time_step << ','
            << next_time_step << ',' << problem.committed_load_factor() << ','
            << nonlinear_iterations;
    const TransientConservationSummary& conservation =
        problem.last_conservation_summary();
    _stream << ',' << conservation.generated_heat_rate << ','
            << conservation.stored_heat_rate << ','
            << conservation.convection_heat_rate << ','
            << conservation.interface_heat_imbalance << ','
            << conservation.dirichlet_heat_input_rate << ','
            << conservation.global_thermal_balance << ','
            << conservation.relative_thermal_balance << ','
            << conservation.unconstrained_thermal_residual_l2 << ','
            << conservation.internal_mechanical_work_increment << ','
            << conservation.pressure_traction_work_increment << ','
            << conservation.dirichlet_reaction_work_increment << ','
            << conservation.contact_work_increment << ','
            << conservation.mechanical_work_balance << ','
            << conservation.relative_mechanical_work_balance << ','
            << conservation.unconstrained_mechanical_residual_l2 << ','
            << conservation.elastic_energy_change << ','
            << conservation.plastic_dissipation_increment << ','
            << conservation.creep_dissipation_increment;
    const std::vector<double>& state = problem.committed_solution();
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        double maximum_temperature =
            -std::numeric_limits<double>::infinity();
        const std::size_t offset = problem.region_node_offset(region);
        for (std::size_t local = 0;
             local < problem.region_mesh(region).nodes().size(); ++local)
            maximum_temperature = std::max(
                maximum_temperature,
                state.at(problem.dof_map().temperature(offset + local)));
        const RegionInelasticSummary history =
            problem.summarize_region_history(region);
        _stream << ',' << maximum_temperature << ','
                << history.maximum_equivalent_plastic_strain << ','
                << history.maximum_equivalent_creep_strain;
    }
    for (std::size_t contact = 0;
         contact < problem.definition().spatial.contacts.size(); ++contact) {
        const InterfaceSummary summary =
            problem.summarize_interface(contact, state);
        _stream << ',' << summary.minimum_gap << ','
                << summary.maximum_contact_pressure << ','
                << summary.total_heat_rate << ','
                << summary.total_contact_force << ','
                << summary.total_tangential_force;
    }
    _stream << '\n';
    _stream.flush();
}

void ExodusResultsIo::write_steady(const std::string& path,
                                   const UnstructuredQuad4Mesh& mesh,
                                   const SteadyProblem& problem,
                                   const std::vector<double>& state) {
    if (path.empty())
        throw std::invalid_argument("Exodus result path must not be empty");
    const std::vector<std::string> nodal =
        nodal_variable_names(problem.definition().contacts);
    const std::vector<std::string> element = stress_variable_names();
    const std::vector<std::string> global =
        global_variable_names(problem.definition().contacts);
    define_variables(path, mesh, nodal, element, global);
    std::vector<std::vector<double>> nodal_values;
    fill_steady_nodal(mesh, problem, state, nodal_values);
    write_step(path, mesh, 1, 1.0, nodal_values,
               steady_elements(mesh, problem, state),
               steady_globals(problem, state));
}

ExodusTransientResultsWriter::ExodusTransientResultsWriter(
    std::string path, UnstructuredQuad4Mesh mesh,
    const TransientProblem& problem)
    : _path(std::move(path)), _mesh(std::move(mesh)),
      _problem_signature(problem.committed_state_signature()), _step_count(0) {
    if (_path.empty())
        throw std::invalid_argument("Exodus result path must not be empty");
    define_variables(
        _path, _mesh,
        nodal_variable_names(problem.definition().spatial.contacts),
        transient_element_variable_names(),
        global_variable_names(problem.definition().spatial.contacts));
}

void ExodusTransientResultsWriter::append(const TransientProblem& problem) {
    if (problem.time_step_active())
        throw std::logic_error(
            "Exodus results cannot be written during an active time step");
    if (problem.committed_state_signature() != _problem_signature)
        throw std::invalid_argument(
            "Exodus result problem does not match writer model");
    std::vector<std::vector<double>> nodal_values;
    fill_transient_nodal(_mesh, problem, nodal_values);
    ++_step_count;
    write_step(_path, _mesh, _step_count, problem.committed_time(),
               nodal_values, transient_elements(_mesh, problem),
               transient_globals(problem));
}

std::size_t ExodusTransientResultsWriter::step_count() const noexcept {
    return _step_count;
}

} // namespace fuelsim
