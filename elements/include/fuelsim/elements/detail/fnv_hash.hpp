#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace fuelsim::detail {
inline constexpr std::uint64_t fnv1a_offset = 14695981039346656037ULL;
inline constexpr std::uint64_t fnv1a_prime = 1099511628211ULL;

inline void fnv1a_bytes(std::uint64_t& hash, const void* data, std::size_t size) noexcept {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= static_cast<std::uint64_t>(bytes[index]);
        hash *= fnv1a_prime;
    }
}

inline std::uint64_t encode_double_bits(double value) noexcept {
    std::uint64_t encoded = 0;
    static_assert(sizeof(encoded) == sizeof(value));
    std::memcpy(&encoded, &value, sizeof(value));
    return encoded;
}

inline double decode_double_bits(std::uint64_t encoded) noexcept {
    double result = 0.0;
    static_assert(sizeof(encoded) == sizeof(result));
    std::memcpy(&result, &encoded, sizeof(result));
    return result;
}
} // namespace fuelsim::detail
