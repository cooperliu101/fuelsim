#pragma once
#include "support/exodus_result_reader.hpp"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace fuelsim::test {
// Plain output records only. These do not reconstruct a Fuelsim problem.
struct OutputPoint final {
    double x, y, z;
};

struct OutputContact final {
    double x, y, z, pressure;
    bool projected, sliding;
    std::array<double, 3> normal_contact_force, tangential_contact_force, tangential_slip;
};

struct OutputContactStep final {
    ExodusResults output;
    std::vector<std::size_t> source_nodes;
    std::vector<OutputContact> contact;
};

inline std::vector<OutputContactStep> read_seven_contact_steps(const std::string& path, std::size_t contact_count) {
    const auto history = read_exodus_nodal_history(path);
    const auto& final = history.back();
    if (final.step_count != 8 || std::abs(final.time - 7.0) > 1e-12)
        throw std::invalid_argument("Friction output requires seven integer-second stages and its initial frame");
    std::vector<OutputContactStep> result;
    for (std::size_t step = 1; step <= 7; ++step) {
        OutputContactStep record;
        record.output = history.at(step);
        const auto& output = record.output;
        if (std::abs(output.time - static_cast<double>(step)) > 1e-12)
            throw std::invalid_argument("Friction path frame time differs from reference");
        const auto& projected = output.nodal("contact_projected_interface");
        for (std::size_t source = 0; source < projected.size(); ++source) {
            if (!std::isfinite(projected[source]))
                continue;
            const auto scalar = [&](const std::string& field) {
                const double value = output.nodal("contact_" + field + "_interface").at(source);
                if (!std::isfinite(value))
                    throw std::invalid_argument("Contact history output is not finite");
                return value;
            };
            record.source_nodes.push_back(source);
            const auto& point = output.nodes.at(source);
            record.contact.push_back({point[0],
                point[1],
                point[2],
                scalar("pressure"),
                scalar("projected") == 1.0,
                scalar("sliding") == 1.0,
                {scalar("normal_force_x"), scalar("normal_force_y"), scalar("normal_force_z")},
                {scalar("tangential_force_x"), scalar("tangential_force_y"), scalar("tangential_force_z")},
                {scalar("total_slip_x"), scalar("total_slip_y"), scalar("total_slip_z")}});
        }
        if (record.contact.size() != contact_count)
            throw std::invalid_argument("Friction path contact-constraint count differs from reference");
        result.push_back(std::move(record));
    }
    return result;
}
} // namespace fuelsim::test
