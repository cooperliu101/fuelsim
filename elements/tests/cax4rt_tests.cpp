#include "support/cax4_common_tests.hpp"

namespace {
using namespace fuelsim::test::cax4;

int run_cax4rt_tests() {
    std::cout << std::scientific << std::setprecision(12);
    return test_cax_kinematics_and_jacobian(true) ? 0 : 1;
}
} // namespace

int main() {
    return run_cax4rt_tests();
}
