#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t step_count = 4;
constexpr std::array<double, 2> friction = {0.3, 0.5};
constexpr double elastic_slip = 0.02;

struct ContactReference final {
    std::size_t step, id, state;
    std::array<double, 3> point, normal_force, tangential_force;
    double slip_1, slip_2, gap, pressure;
};

struct ReactionReference final {
    std::size_t step;
    std::array<std::array<double, 3>, 2> pair;
    std::array<double, 3> global;
};

struct OutputContact final {
    std::array<double, 3> normal_contact_force, tangential_contact_force, displacement;
    double gap, pressure, contact_force, primary_face;
    bool sliding, projected;
};

struct OutputHistory final {
    std::array<double, 3> cartesian_elastic_tangential_slip;
};

struct StepState final {
    std::array<std::vector<OutputContact>, 2> contact;
    std::array<std::vector<OutputHistory>, 2> histories;
};

bool check(bool condition, const std::string& message) {
    if (!condition)
        std::cerr << "[FAIL] " << message << '\n';
    return condition;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ','))
        result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    if (column >= values.size())
        throw std::invalid_argument("Incomplete B4.8 CSV row: " + input_path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid B4.8 CSV number: " + input_path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& input_path) {
    const double value = number(values, column, input_path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max())
        || std::floor(value) != value)
        throw std::invalid_argument("Invalid B4.8 CSV index: " + input_path);
    return static_cast<std::size_t>(value);
}

std::vector<ContactReference> read_contact(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input)
        throw std::runtime_error("Could not read B4.8 contact reference: " + input_path);
    std::string line;
    if (!std::getline(input, line)
        || line
               != "step,id,x,y,z,normal_x,normal_y,normal_z,tangential_x,tangential_y,"
                  "tangential_z,slip1,slip2,gap,pressure,state")
        throw std::invalid_argument("Unexpected B4.8 contact CSV header: " + input_path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split(line);
        result.push_back({index_value(values, 0, input_path),
            index_value(values, 1, input_path),
            index_value(values, 15, input_path),
            {number(values, 2, input_path), number(values, 3, input_path), number(values, 4, input_path)},
            {number(values, 5, input_path), number(values, 6, input_path), number(values, 7, input_path)},
            {number(values, 8, input_path), number(values, 9, input_path), number(values, 10, input_path)},
            number(values, 11, input_path),
            number(values, 12, input_path),
            number(values, 13, input_path),
            number(values, 14, input_path)});
    }
    if (result.size() != step_count * 4)
        throw std::invalid_argument("B4.8 contact reference must have sixteen rows");
    return result;
}

std::vector<ReactionReference> read_reactions(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input)
        throw std::runtime_error("Could not read B4.8 reaction reference: " + input_path);
    std::string line;
    if (!std::getline(input, line)
        || line != "step,pair_a_x,pair_a_y,pair_a_z,pair_b_x,pair_b_y,pair_b_z,global_x,global_y,global_z")
        throw std::invalid_argument("Unexpected B4.8 reaction CSV header: " + input_path);
    std::vector<ReactionReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split(line);
        result.push_back({index_value(values, 0, input_path),
            {{{number(values, 1, input_path), number(values, 2, input_path), number(values, 3, input_path)},
                {number(values, 4, input_path), number(values, 5, input_path), number(values, 6, input_path)}}},
            {number(values, 7, input_path), number(values, 8, input_path), number(values, 9, input_path)}});
    }
    if (result.size() != step_count)
        throw std::invalid_argument("B4.8 reaction reference must have four rows");
    return result;
}

std::vector<double> read_energy(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input)
        throw std::runtime_error("Could not read B4.8 energy reference: " + input_path);
    std::string line;
    if (!std::getline(input, line) || line != "step,allfd")
        throw std::invalid_argument("Unexpected B4.8 energy CSV header: " + input_path);
    std::vector<double> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split(line);
        if (index_value(values, 0, input_path) != result.size() + 1)
            throw std::invalid_argument("B4.8 energy steps are not consecutive");
        result.push_back(number(values, 1, input_path));
    }
    if (result.size() != step_count)
        throw std::invalid_argument("B4.8 energy reference must have four rows");
    return result;
}

std::array<double, 3> reference_slip(const ContactReference& value) {
    return {0.0, -value.slip_2, value.slip_1};
}

std::array<double, 3> reference_elastic_slip(const ContactReference& value, double mu) {
    std::array<double, 3> result{};
    const double normal =
        std::sqrt(value.normal_force[0] * value.normal_force[0] + value.normal_force[1] * value.normal_force[1]
                  + value.normal_force[2] * value.normal_force[2]);
    if (normal == 0.0 || mu == 0.0)
        return result;
    for (std::size_t component = 0; component < 3; ++component)
        result[component] = -value.tangential_force[component] * elastic_slip / (mu * normal);
    return result;
}

bool compare(const std::vector<std::array<double, 3>>& mesh,
    const std::array<std::vector<std::size_t>, 2>& sources,
    const std::vector<StepState>& states,
    const std::array<std::vector<ContactReference>, 2>& references,
    const std::vector<ReactionReference>& reactions,
    const std::vector<double>& allfd) {
    if (states.size() != step_count)
        throw std::invalid_argument("B4.8 Fuelsim path must contain four states");
    std::array<fuelsim::test::FieldErrorMetrics, 2> normal, tangent_y, tangent_z, slip_y, slip_z, gap, pressure,
        resultant, force_center, dissipation;
    fuelsim::test::FieldErrorMetrics global_reaction, global_dissipation;
    std::array<std::vector<std::array<double, 3>>, 2> previous_actual_plastic, previous_reference_plastic;
    for (std::size_t pair = 0; pair < 2; ++pair) {
        previous_actual_plastic[pair].resize(sources[pair].size());
        previous_reference_plastic[pair].resize(sources[pair].size());
    }
    std::array<double, 2> actual_cumulative_dissipation{}, reference_cumulative_dissipation{};
    double maximum_coordinate_difference = 0.0, maximum_zero_reference_difference = 0.0,
           maximum_action_reaction_difference = 0.0, maximum_abaqus_shear_normal_component = 0.0,
           maximum_abaqus_shear_normal_resultant = 0.0;
    bool all_active = true, first_sticks = true, later_slides = true, crossed_different_primary_faces = false,
         states_match = true;
    std::array<double, 2> first_primary_face = {-1, -1};
    for (std::size_t step = 0; step < states.size(); ++step) {
        std::array<std::array<double, 3>, 2> actual_resultant{};
        std::array<double, 3> actual_global{};
        for (std::size_t pair = 0; pair < 2; ++pair) {
            std::array<double, 3> actual_center_sum{}, reference_center_sum{};
            double actual_center_weight = 0.0, reference_center_weight = 0.0;
            double actual_increment_dissipation = 0.0, reference_increment_dissipation = 0.0;
            double abaqus_shear_normal_resultant = 0.0;
            std::size_t sliding_nodes = 0;
            for (std::size_t node = 0; node < states[step].contact[pair].size(); ++node) {
                const auto found =
                    std::find_if(references[pair].begin(), references[pair].end(), [&](const ContactReference& value) {
                        return value.step == step + 1 && value.id == sources[pair][node];
                    });
                if (found == references[pair].end())
                    throw std::invalid_argument("B4.8 contact mapping is incomplete");
                const auto& actual = states[step].contact[pair][node];
                const auto& point = mesh.at(sources[pair][node]);
                maximum_coordinate_difference = std::max({maximum_coordinate_difference,
                    std::abs(point[0] - found->point[0]),
                    std::abs(point[1] - found->point[1]),
                    std::abs(point[2] - found->point[2])});
                normal[pair].add(-actual.normal_contact_force[0], found->normal_force[0]);
                tangent_y[pair].add(-actual.tangential_contact_force[1], found->tangential_force[1]);
                tangent_z[pair].add(-actual.tangential_contact_force[2], found->tangential_force[2]);
                const auto& displacement = actual.displacement;
                slip_y[pair].add(displacement[1], -found->slip_2);
                slip_z[pair].add(displacement[2], found->slip_1);
                gap[pair].add(actual.gap, found->gap);
                pressure[pair].add(actual.pressure, found->pressure);
                maximum_zero_reference_difference = std::max({maximum_zero_reference_difference,
                    std::abs(actual.normal_contact_force[1]),
                    std::abs(actual.normal_contact_force[2]),
                    std::abs(actual.tangential_contact_force[0]),
                    std::abs(found->normal_force[1]),
                    std::abs(found->normal_force[2])});
                maximum_abaqus_shear_normal_component =
                    std::max(maximum_abaqus_shear_normal_component, std::abs(found->tangential_force[0]));
                abaqus_shear_normal_resultant += found->tangential_force[0];
                const std::array<double, 3> actual_force = {-actual.normal_contact_force[0]
                                                                - actual.tangential_contact_force[0],
                    -actual.normal_contact_force[1] - actual.tangential_contact_force[1],
                    -actual.normal_contact_force[2] - actual.tangential_contact_force[2]};
                const std::array<double, 3> actual_current = {point[0] + displacement[0],
                    point[1] + displacement[1],
                    point[2] + displacement[2]};
                const std::array<double, 3> reference_current = actual_current;
                for (std::size_t component = 0; component < 3; ++component) {
                    actual_resultant[pair][component] += actual_force[component];
                    actual_global[component] += actual_force[component];
                    actual_center_sum[component] += actual.contact_force * actual_current[component];
                    reference_center_sum[component] += std::abs(found->normal_force[0]) * reference_current[component];
                }
                actual_center_weight += actual.contact_force;
                reference_center_weight += std::abs(found->normal_force[0]);
                const std::array<double, 3> reference_total_slip = reference_slip(*found),
                                            reference_elastic = reference_elastic_slip(*found, friction[pair]);
                std::array<double, 3> actual_plastic{}, reference_plastic{};
                for (std::size_t component = 0; component < 3; ++component) {
                    const double actual_total_slip = component == 0 ? 0.0 : displacement[component];
                    actual_plastic[component] =
                        actual_total_slip
                        - states[step].histories[pair][node].cartesian_elastic_tangential_slip[component];
                    reference_plastic[component] = reference_total_slip[component] - reference_elastic[component];
                    if (component > 0) {
                        actual_increment_dissipation -=
                            (-actual.tangential_contact_force[component])
                            * (actual_plastic[component] - previous_actual_plastic[pair][node][component]);
                        reference_increment_dissipation -=
                            found->tangential_force[component]
                            * (reference_plastic[component] - previous_reference_plastic[pair][node][component]);
                    }
                }
                previous_actual_plastic[pair][node] = actual_plastic;
                previous_reference_plastic[pair][node] = reference_plastic;
                const std::size_t actual_state = actual.pressure <= 0.0 ? 0 : (actual.sliding ? 2 : 1);
                states_match = states_match && actual_state == found->state;
                all_active = all_active && actual.projected && actual.pressure > 0.0 && found->pressure > 0.0;
                first_sticks = first_sticks && (step != 0 || (!actual.sliding && found->state == 1));
                if (actual.sliding)
                    ++sliding_nodes;
                if (step == 0 && node == 0)
                    first_primary_face[pair] = actual.primary_face;
                if (step > 0 && actual.primary_face != first_primary_face[pair])
                    crossed_different_primary_faces = true;
            }
            later_slides = later_slides && (step == 0 || sliding_nodes > 0);
            maximum_abaqus_shear_normal_resultant =
                std::max(maximum_abaqus_shear_normal_resultant, std::abs(abaqus_shear_normal_resultant));
            actual_cumulative_dissipation[pair] += actual_increment_dissipation;
            reference_cumulative_dissipation[pair] += reference_increment_dissipation;
            dissipation[pair].add(actual_cumulative_dissipation[pair],
                step == 0 ? 0.0 : reference_cumulative_dissipation[pair]);
            for (std::size_t component = 0; component < 3; ++component) {
                resultant[pair].add(actual_resultant[pair][component], -reactions[step].pair[pair][component]);
                force_center[pair].add(actual_center_sum[component] / actual_center_weight,
                    reference_center_sum[component] / reference_center_weight);
                maximum_action_reaction_difference = std::max(maximum_action_reaction_difference,
                    std::abs(actual_resultant[pair][component] + reactions[step].pair[pair][component]));
            }
        }
        const double actual_total_dissipation = actual_cumulative_dissipation[0] + actual_cumulative_dissipation[1];
        global_dissipation.add(actual_total_dissipation, allfd[step]);
        for (std::size_t component = 0; component < 3; ++component)
            global_reaction.add(actual_global[component], -reactions[step].global[component]);
    }
    const auto print = [](const std::string& name, const fuelsim::test::FieldErrorMetrics& metric) {
        if (metric.has_relative_norm())
            fuelsim::test::print_relative_metrics(name, metric);
        else
            fuelsim::test::print_absolute_metrics(name, metric);
    };
    for (std::size_t pair = 0; pair < 2; ++pair) {
        const std::string prefix = std::string("b48_pair_") + static_cast<char>('a' + pair) + '_';
        print(prefix + "normal_force_x", normal[pair]);
        print(prefix + "tangential_force_y", tangent_y[pair]);
        print(prefix + "tangential_force_z", tangent_z[pair]);
        print(prefix + "slip_y", slip_y[pair]);
        print(prefix + "slip_z", slip_z[pair]);
        print(prefix + "gap", gap[pair]);
        print(prefix + "pressure", pressure[pair]);
        print(prefix + "resultant", resultant[pair]);
        print(prefix + "force_center", force_center[pair]);
        print(prefix + "friction_dissipation", dissipation[pair]);
    }
    print("b48_global_reaction", global_reaction);
    print("b48_global_friction_dissipation", global_dissipation);
    std::cout << "b48_maximum_mesh_coordinate_difference=" << maximum_coordinate_difference << '\n'
              << "b48_maximum_zero_reference_force_difference=" << maximum_zero_reference_difference << '\n'
              << "b48_maximum_action_reaction_difference=" << maximum_action_reaction_difference << '\n'
              << "b48_abaqus_maximum_nodal_shear_normal_component=" << maximum_abaqus_shear_normal_component << '\n'
              << "b48_abaqus_maximum_shear_normal_resultant=" << maximum_abaqus_shear_normal_resultant << '\n';
    constexpr double tolerance = 1.0e-2, zero_tolerance = 1.0e-7;
    const auto passes = [&](const fuelsim::test::FieldErrorMetrics& metric) {
        return (!metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, tolerance))
               && metric.maximum_zero_reference_difference < zero_tolerance;
    };
    bool fields_pass = passes(global_reaction) && passes(global_dissipation);
    for (std::size_t pair = 0; pair < 2; ++pair)
        fields_pass = fields_pass && passes(normal[pair]) && passes(tangent_y[pair]) && passes(tangent_z[pair])
                      && passes(slip_y[pair]) && passes(slip_z[pair]) && passes(gap[pair]) && passes(pressure[pair])
                      && passes(resultant[pair]) && passes(force_center[pair]) && passes(dissipation[pair]);
    return check(maximum_coordinate_difference < 3.0e-8, "B4.8 Abaqus and Fuelsim use the same tracked mesh")
           && check(all_active && first_sticks && later_slides && crossed_different_primary_faces,
               "B4.8 keeps both pairs active, starts in sticking, slides both pairs, and transfers primary-face "
               "ownership")
           && check(states_match,
               "B4.8 open, sticking, and sliding states agree with Abaqus at every accepted increment")
           && check(maximum_zero_reference_difference < zero_tolerance,
               "B4.8 theoretical-zero contact-force components pass the separate absolute check")
           && check(maximum_abaqus_shear_normal_resultant < 1.0e-10,
               "B4.8 Abaqus nodal shear leakage normal to the plane cancels in each pair resultant")
           && check(fields_pass,
               "B4.8 all per-pair nodal fields, resultants, force centers, friction dissipation, and "
               "global reactions pass the one-percent metrics");
}

} // namespace

namespace fuelsim::test {
bool check_hex8_multi_contact_path(const std::string& output, const std::string& reference_directory) {
    const auto frames = read_exodus_nodal_history(output);
    if (frames.size() != 5)
        throw std::invalid_argument("B4.8 requires initial state and four increments");
    const std::array<std::vector<ContactReference>, 2> refs = {
        read_contact(reference_directory + "/b48_hex8_multi_contact_pair_a.csv"),
        read_contact(reference_directory + "/b48_hex8_multi_contact_pair_b.csv")};
    const auto reactions = read_reactions(reference_directory + "/b48_hex8_multi_contact_reaction.csv");
    const auto energy = read_energy(reference_directory + "/b48_hex8_multi_contact_energy.csv");
    std::array<std::vector<std::size_t>, 2> sources;
    for (std::size_t pair = 0; pair < 2; ++pair) {
        std::set<std::pair<std::size_t, std::size_t>> seen;
        for (const auto& ref : refs[pair]) {
            if (ref.step < 1 || ref.step > 4 || ref.id >= frames[0].nodes.size()
                || !seen.emplace(ref.step, ref.id).second)
                throw std::invalid_argument("B4.8 reference association is invalid");
            if (ref.step == 1)
                sources[pair].push_back(ref.id);
        }
        std::sort(sources[pair].begin(), sources[pair].end());
        if (sources[pair].size() != 4)
            throw std::invalid_argument("B4.8 must contain four constraints per pair");
    }
    std::vector<StepState> states;
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        if (frame.time != static_cast<double>(step) || reactions[step - 1].step != step)
            throw std::invalid_argument("B4.8 increment association is invalid");
        StepState state;
        for (std::size_t pair = 0; pair < 2; ++pair) {
            const std::string suffix = pair == 0 ? "_pair_a" : "_pair_b";
            const auto& flags = frame.nodal("contact_projected" + suffix);
            if (std::count_if(flags.begin(), flags.end(), [](double value) { return !std::isnan(value); }) != 4)
                throw std::invalid_argument("B4.8 output must cover exactly four constraints per pair");
            for (auto n : sources[pair]) {
                OutputContact contact{};
                OutputHistory history{};
                contact.gap = frame.nodal("contact_gap" + suffix).at(n);
                contact.pressure = frame.nodal("contact_pressure" + suffix).at(n);
                contact.primary_face = frame.nodal("contact_primary_face" + suffix).at(n);
                contact.sliding = frame.nodal("contact_sliding" + suffix).at(n) == 1;
                contact.projected = flags.at(n) == 1;
                for (std::size_t c = 0; c < 3; ++c) {
                    const std::string axis(1, "xyz"[c]);
                    contact.normal_contact_force[c] = frame.nodal("contact_normal_force_" + axis + suffix).at(n);
                    contact.tangential_contact_force[c] =
                        frame.nodal("contact_tangential_force_" + axis + suffix).at(n);
                    contact.displacement[c] = frame.nodal("displacement_" + axis).at(n);
                    history.cartesian_elastic_tangential_slip[c] =
                        frame.nodal("contact_elastic_slip_" + axis + suffix).at(n);
                }
                contact.contact_force = std::hypot(contact.normal_contact_force[0],
                    contact.normal_contact_force[1],
                    contact.normal_contact_force[2]);
                state.contact[pair].push_back(contact);
                state.histories[pair].push_back(history);
            }
        }
        states.push_back(state);
    }
    return compare(frames[0].nodes, sources, states, refs, reactions, energy);
}
} // namespace fuelsim::test
