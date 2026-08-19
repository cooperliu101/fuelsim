#include "fuelsim/io/checkpoint.hpp"
#include "core/problem_backend_access.hpp"
#include "detail/fnv_hash.hpp"
#include "fuelsim/core/transient_problem.hpp"
#include "io/problem_signature.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
namespace {
constexpr std::array<unsigned char, 16> checkpoint_magic = {
    'F', 'U', 'E', 'L', 'S', 'I', 'M', '_', 'C', 'H', 'E', 'C', 'K', 'P', 'T', '\0'};
constexpr std::uint32_t checkpoint_version = 11U;
constexpr std::uint32_t endian_marker = 0x01020304U;
constexpr std::uint64_t maximum_checkpoint_bytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;

std::uint64_t checksum(const std::vector<unsigned char>& bytes) {
    std::uint64_t value = detail::fnv1a_offset;
    detail::fnv1a_bytes(value, bytes.data(), bytes.size());
    return value;
}

class BinaryBuffer final {
  public:
    void append_u32(std::uint32_t value) {
        for (std::size_t byte = 0; byte < 4; ++byte) _bytes.push_back(static_cast<unsigned char>(value >> (8U * byte)));
    }

    void append_u64(std::uint64_t value) {
        for (std::size_t byte = 0; byte < 8; ++byte) _bytes.push_back(static_cast<unsigned char>(value >> (8U * byte)));
    }

    void append_double(double value) { append_u64(detail::encode_double_bits(value)); }

    void append_bytes(const unsigned char* data, std::size_t size) { _bytes.insert(_bytes.end(), data, data + size); }

    const std::vector<unsigned char>& bytes() const noexcept { return _bytes; }

  private:
    std::vector<unsigned char> _bytes;
};

class BinaryCursor final {
  public:
    explicit BinaryCursor(const std::vector<unsigned char>& bytes) : _bytes(bytes), _position(0) {}

    std::uint32_t read_u32() {
        require(4);
        std::uint32_t result = 0;
        for (std::size_t byte = 0; byte < 4; ++byte)
            result |= static_cast<std::uint32_t>(_bytes[_position++]) << (8U * byte);
        return result;
    }

    std::uint64_t read_u64() {
        require(8);
        std::uint64_t result = 0;
        for (std::size_t byte = 0; byte < 8; ++byte)
            result |= static_cast<std::uint64_t>(_bytes[_position++]) << (8U * byte);
        return result;
    }

    double read_double() { return detail::decode_double_bits(read_u64()); }

    void read_bytes(unsigned char* destination, std::size_t size) {
        require(size);
        std::memcpy(destination, _bytes.data() + _position, size);
        _position += size;
    }

    bool at_end() const noexcept { return _position == _bytes.size(); }

  private:
    void require(std::size_t size) const {
        if (size > _bytes.size() - _position) throw std::runtime_error("Checkpoint payload is truncated");
    }

    const std::vector<unsigned char>& _bytes;
    std::size_t _position;
};

void append_material_values(BinaryBuffer& payload, const double* elastic, const double* plastic, const double* creep,
    std::size_t count, double equivalent_plastic, double equivalent_creep, const double* stress) {
    for (const double* values : {elastic, plastic, creep})
        for (std::size_t component = 0; component < count; ++component) payload.append_double(values[component]);
    payload.append_double(equivalent_plastic);
    payload.append_double(equivalent_creep);
    for (std::size_t component = 0; component < count; ++component) payload.append_double(stress[component]);
}

void read_material_values(BinaryCursor& payload, double* elastic, double* plastic, double* creep, std::size_t count,
    double& equivalent_plastic, double& equivalent_creep, double* stress) {
    for (double* values : {elastic, plastic, creep})
        for (std::size_t component = 0; component < count; ++component) values[component] = payload.read_double();
    equivalent_plastic = payload.read_double();
    equivalent_creep = payload.read_double();
    for (std::size_t component = 0; component < count; ++component) stress[component] = payload.read_double();
}

void append_material_point(BinaryBuffer& payload, const MaterialPointState& state) {
    const double stress[] = {state.stress.rr, state.stress.zz, state.stress.hoop, state.stress.rz};
    append_material_values(payload, state.elastic_strain.data(), state.plastic_strain.data(), state.creep_strain.data(),
        state.elastic_strain.size(), state.equivalent_plastic_strain, state.equivalent_creep_strain, stress);
}

void read_material_point(BinaryCursor& payload, MaterialPointState& state) {
    double stress[4];
    read_material_values(payload, state.elastic_strain.data(), state.plastic_strain.data(), state.creep_strain.data(),
        4, state.equivalent_plastic_strain, state.equivalent_creep_strain, stress);
    state.stress = {stress[0], stress[1], stress[2], stress[3]};
}

void append_material_point(BinaryBuffer& payload, const CartesianMaterialPointState& state) {
    const double stress[] = {
        state.stress.xx, state.stress.yy, state.stress.zz, state.stress.xy, state.stress.yz, state.stress.xz};
    append_material_values(payload, state.elastic_strain.data(), state.plastic_strain.data(), state.creep_strain.data(),
        state.elastic_strain.size(), state.equivalent_plastic_strain, state.equivalent_creep_strain, stress);
}

void read_material_point(BinaryCursor& payload, CartesianMaterialPointState& state) {
    double stress[6];
    read_material_values(payload, state.elastic_strain.data(), state.plastic_strain.data(), state.creep_strain.data(),
        6, state.equivalent_plastic_strain, state.equivalent_creep_strain, stress);
    state.stress = {stress[0], stress[1], stress[2], stress[3], stress[4], stress[5]};
}

void append_conservation(BinaryBuffer& payload, const TransientConservationSummary& summary) {
    for (const TransientConservationField& field : transient_conservation_fields)
        payload.append_double(summary.*field.member);
}

TransientConservationSummary read_conservation(BinaryCursor& payload) {
    TransientConservationSummary result;
    for (const TransientConservationField& field : transient_conservation_fields) {
        result.*field.member = payload.read_double();
        if (!std::isfinite(result.*field.member))
            throw std::runtime_error("Checkpoint conservation summary is invalid");
    }
    return result;
}

BinaryBuffer state_payload(const TransientProblem& problem, double next_time_step) {
    BinaryBuffer payload;
    payload.append_u64(transient_problem_signature(problem));
    const bool cartesian = problem.is_cartesian_3d();
    const TransientCommittedState state = BackendAccess::committed_state(problem);
    payload.append_double(state.time);
    payload.append_double(state.load_factor);
    payload.append_double(next_time_step);
    append_conservation(payload, state.conservation);
    for (const double value : state.solution) payload.append_double(value);
    for (const auto& contact : state.contact_histories) {
        for (const ContactPointHistory& history : contact) {
            payload.append_double(history.elastic_tangential_slip);
            payload.append_u32(history.sliding ? 1U : 0U);
            payload.append_double(history.normal_multiplier);
            for (double component : history.cartesian_elastic_tangential_slip) payload.append_double(component);
        }
    }
    if (cartesian) {
        for (const auto& region : state.cartesian_material_histories)
            for (const Hex8MaterialHistory& element : region)
                for (const CartesianMaterialPointState& point : element) append_material_point(payload, point);
        return payload;
    }
    for (std::size_t region = 0; region < state.material_histories.size(); ++region) {
        for (std::size_t element = 0; element < state.material_histories[region].size(); ++element)
            for (std::size_t q = 0; q < 4; ++q)
                append_material_point(payload, state.material_histories[region][element][q]);
    }
    return payload;
}

void append_checkpoint_header(BinaryBuffer& file, const std::vector<unsigned char>& payload) {
    file.append_bytes(checkpoint_magic.data(), checkpoint_magic.size());
    file.append_u32(checkpoint_version);
    file.append_u32(endian_marker);
    file.append_u64(static_cast<std::uint64_t>(payload.size()));
    file.append_u64(checksum(payload));
}

std::vector<unsigned char> read_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("Could not open checkpoint '" + path + "'");
    const std::streamoff length = input.tellg();
    if (length < 0 || static_cast<std::uint64_t>(length) > maximum_checkpoint_bytes)
        throw std::runtime_error("Checkpoint size is invalid: " + path);
    input.seekg(0, std::ios::beg);
    std::vector<unsigned char> result(static_cast<std::size_t>(length));
    if (!result.empty()) input.read(reinterpret_cast<char*>(result.data()), length);
    if (!input) throw std::runtime_error("Could not read checkpoint '" + path + "'");
    return result;
}

void write_atomic(const std::string& path, const std::vector<unsigned char>& bytes) {
    const std::string temporary = path + ".tmp";
    try {
        std::ofstream output(temporary, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!output) throw std::runtime_error("Could not create checkpoint '" + temporary + "'");
        output.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        output.close();
        if (!output) throw std::runtime_error("Could not finish checkpoint '" + temporary + "'");
        if (std::rename(temporary.c_str(), path.c_str()) != 0)
            throw std::runtime_error("Could not install checkpoint '" + path + "': " + std::strerror(errno));
    } catch (...) {
        (void)std::remove(temporary.c_str());
        throw;
    }
}
} // namespace

void write_transient_checkpoint(const std::string& path, const TransientProblem& problem, double next_time_step) {
    if (path.empty()) throw std::invalid_argument("Checkpoint path must not be empty");
    if (problem.time_step_active()) throw std::logic_error("Checkpoint cannot be written during an active time step");
    if (!std::isfinite(next_time_step) || !(next_time_step > 0.0))
        throw std::invalid_argument("Checkpoint next time step must be finite and positive");
    const BinaryBuffer payload = state_payload(problem, next_time_step);
    BinaryBuffer file;
    append_checkpoint_header(file, payload.bytes());
    file.append_bytes(payload.bytes().data(), payload.bytes().size());
    write_atomic(path, file.bytes());
}

double restore_transient_checkpoint(const std::string& path, TransientProblem& problem) {
    const std::vector<unsigned char> file = read_file(path);
    BinaryCursor header(file);
    std::array<unsigned char, checkpoint_magic.size()> magic{};
    header.read_bytes(magic.data(), magic.size());
    if (magic != checkpoint_magic) throw std::runtime_error("Checkpoint magic does not match fuelsim");
    if (header.read_u32() != checkpoint_version) throw std::runtime_error("Checkpoint version is not supported");
    if (header.read_u32() != endian_marker) throw std::runtime_error("Checkpoint endian marker is invalid");
    const std::uint64_t payload_size = header.read_u64();
    const std::uint64_t expected_checksum = header.read_u64();
    constexpr std::size_t header_size = 16U + 4U + 4U + 8U + 8U;
    if (payload_size != file.size() - header_size) throw std::runtime_error("Checkpoint payload length is invalid");
    std::vector<unsigned char> payload_bytes(file.begin() + static_cast<std::ptrdiff_t>(header_size), file.end());
    if (checksum(payload_bytes) != expected_checksum)
        throw std::runtime_error("Checkpoint payload checksum does not match");
    BinaryCursor payload(payload_bytes);
    if (payload.read_u64() != transient_problem_signature(problem))
        throw std::runtime_error("Checkpoint model signature does not match the current problem");
    const bool cartesian = problem.is_cartesian_3d();
    TransientCommittedState state;
    state.time = payload.read_double();
    state.load_factor = payload.read_double();
    const double next_time_step = payload.read_double();
    state.conservation = read_conservation(payload);
    if (!std::isfinite(next_time_step) || !(next_time_step > 0.0))
        throw std::runtime_error("Checkpoint next time step is invalid");
    state.solution.resize(problem.dof_count());
    for (double& value : state.solution) value = payload.read_double();
    const TransientCommittedState expected_state = BackendAccess::committed_state(problem);
    const auto& expected_histories = expected_state.contact_histories;
    state.contact_histories.resize(expected_histories.size());
    for (std::size_t contact = 0; contact < expected_histories.size(); ++contact) {
        state.contact_histories[contact].resize(expected_histories[contact].size());
        for (ContactPointHistory& history : state.contact_histories[contact]) {
            history.elastic_tangential_slip = payload.read_double();
            const std::uint32_t sliding = payload.read_u32();
            if (sliding > 1U) throw std::runtime_error("Checkpoint friction state is invalid");
            history.sliding = sliding == 1U;
            history.normal_multiplier = payload.read_double();
            if (!std::isfinite(history.normal_multiplier) || history.normal_multiplier < 0.0)
                throw std::runtime_error("Checkpoint normal contact multiplier is invalid");
            for (double& component : history.cartesian_elastic_tangential_slip) {
                component = payload.read_double();
                if (!std::isfinite(component))
                    throw std::runtime_error("Checkpoint Cartesian friction state is invalid");
            }
        }
    }
    if (cartesian) {
        const auto& expected_histories = BackendAccess::cartesian_material_histories(problem);
        state.cartesian_material_histories.resize(expected_histories.size());
        for (std::size_t region = 0; region < expected_histories.size(); ++region) {
            state.cartesian_material_histories[region].resize(expected_histories[region].size());
            for (Hex8MaterialHistory& element : state.cartesian_material_histories[region])
                for (CartesianMaterialPointState& point : element) read_material_point(payload, point);
        }
        if (!payload.at_end()) throw std::runtime_error("Checkpoint payload contains trailing data");
        BackendAccess::restore_committed_state(problem, std::move(state));
        return next_time_step;
    }
    const rz::TransientBackendView backend = BackendAccess::transient(problem);
    state.material_histories.resize(backend.spatial.region_count());
    for (std::size_t region = 0; region < backend.spatial.region_count(); ++region) {
        const std::size_t elements = backend.spatial.region_mesh(region).elements().size();
        state.material_histories[region].resize(elements);
        for (std::size_t element = 0; element < elements; ++element)
            for (std::size_t q = 0; q < 4; ++q)
                read_material_point(payload, state.material_histories[region][element][q]);
    }
    if (!payload.at_end()) throw std::runtime_error("Checkpoint payload contains trailing data");
    BackendAccess::restore_committed_state(problem, std::move(state));
    return next_time_step;
}
} // namespace fuelsim
