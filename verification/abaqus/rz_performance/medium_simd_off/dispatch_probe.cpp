#include "derivative_kernels.hpp"
#include <iostream>

int main() {
    namespace k = adlite::detail::derivative_kernels;
    for (std::size_t width : {5u, 6u}) {
        const auto& f = k::functions_for_width(width);
        std::cout << "width=" << width << " scale_scalar=" << (f.scale == k::scale_scalar)
                  << " accumulate_scalar=" << (f.accumulate_scaled == k::accumulate_scaled_scalar)
                  << " combine_scalar=" << (f.combine == k::combine_scalar) << '\n';
    }
}
