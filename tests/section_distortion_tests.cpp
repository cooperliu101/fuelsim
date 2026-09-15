#include "quad8_face.hpp"
#include "solver/petsc_solver.hpp"
#include "solver/section_distortion.hpp"
#include "support/section_fixture.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
std::array<double, 2>
sample(const fuelsim::CrossSection& section, const fuelsim::SectionDistortionMode& mode, double x, double y) {
    for (std::size_t cell = 0; cell < section.points().size(); cell += 9) {
        const auto& point = section.points()[cell];
        const auto& lower = section.nodes()[point.nodes[0]];
        const auto& upper = section.nodes()[point.nodes[2]];
        if (x < lower.x - 1.0e-14 || x > upper.x + 1.0e-14 || y < lower.y - 1.0e-14 || y > upper.y + 1.0e-14)
            continue;
        fuelsim::Quad8FaceCoordinates coordinates{};
        for (std::size_t i = 0; i < 8; ++i)
            coordinates[i] = section.nodes()[point.nodes[i]];
        const auto shape = fuelsim::make_quad8_face_mechanical_point(coordinates,
            2.0 * (x - lower.x) / (upper.x - lower.x) - 1.0,
            2.0 * (y - lower.y) / (upper.y - lower.y) - 1.0,
            1.0)
                               .displacement_shape;
        std::array<double, 2> value{};
        for (std::size_t i = 0; i < 8; ++i) {
            value[0] += shape[i] * mode.transverse[point.nodes[i]];
            value[1] += shape[i] * mode.transverse[section.nodes().size() + point.nodes[i]];
        }
        return value;
    }
    throw std::runtime_error("Mode comparison sample is outside the section");
}

void verify() {
    std::vector<double> previous;
    std::unique_ptr<fuelsim::CrossSection> previous_section;
    fuelsim::SectionDistortionSpectrum previous_spectrum;
    for (int refinement : {1, 2, 3}) {
        std::vector<double> x;
        for (int i = 0; i <= refinement; ++i)
            x.push_back(-0.001 + 0.002 * static_cast<double>(i) / static_cast<double>(refinement));
        const auto section =
            fuelsim::test::rectangle(x, 0.02, {fuelsim::test::region("al", 70e9, 0.3)}, false, 4 * refinement);
        const auto spectrum = fuelsim::build_section_distortion_modes(section, 6, 4);
        for (std::size_t i = 0; i < spectrum.shear_free_modes.size(); ++i) {
            const auto& mode = spectrum.shear_free_modes[i];
            std::cout << "shear_free_mode=" << i << " KD_over_W=" << mode.distortion_energy
                      << " W_norm=" << mode.section_norm << " residual=" << mode.eigen_relative_residual
                      << " orthogonality=" << mode.orthogonality_error << '\n';
            if (!std::isinf(mode.eigenvalue) || std::abs(mode.section_norm - 1.0) > 1.0e-9
                || mode.classic_projection > 1.0e-9 || mode.orthogonality_error > 1.0e-7
                || mode.eigen_relative_residual > 1.0e-7
                || (i != 0 && mode.distortion_energy < spectrum.shear_free_modes[i - 1].distortion_energy))
                throw std::runtime_error("Shear-free normalization, ordering or operator residual failed");
            // Direct quadrature of gamma verifies cancellation independently of
            // the condensed matrix and its transformed eigensystem.
            double shear = 0.0, transverse = 0.0;
            for (const auto& point : section.points()) {
                std::array<double, 2> phi{}, gradient{};
                for (std::size_t node = 0; node < 8; ++node)
                    for (std::size_t c = 0; c < 2; ++c) {
                        phi[c] += point.shape[node] * mode.transverse[c * section.nodes().size() + point.nodes[node]];
                        gradient[c] += point.gradient[node][c] * mode.warping.axial[point.nodes[node]];
                    }
                for (std::size_t c = 0; c < 2; ++c) {
                    shear += point.weight * std::pow(phi[c] + gradient[c], 2);
                    transverse += point.weight * phi[c] * phi[c];
                }
            }
            if (shear / transverse > 1.0e-16)
                throw std::runtime_error("Shear kernel did not cancel pointwise shear");
        }
        std::cout << "distortion_refinement=" << refinement << " nodes=" << section.nodes().size()
                  << " shear_kernel=" << spectrum.shear_kernel_dimension << '\n';
        for (std::size_t i = 0; i < spectrum.modes.size(); ++i) {
            const auto& mode = spectrum.modes[i];
            std::cout << "mode=" << i << " lambda=" << mode.eigenvalue << " norm=" << mode.shear_norm
                      << " W_norm=" << mode.section_norm << " orthogonality=" << mode.orthogonality_error
                      << " classic_projection=" << mode.classic_projection
                      << " residual=" << mode.eigen_relative_residual << '\n';
            if (std::abs(mode.shear_norm - 1.0) > 1.0e-9 || !std::isfinite(mode.eigenvalue)
                || (i > 0 && mode.eigenvalue < spectrum.modes[i - 1].eigenvalue))
                throw std::runtime_error("Distortion normalization or ordering failed");
            if (refinement == 3 && i < 3 && std::abs(mode.eigenvalue / previous[i] - 1.0) > 0.01)
                throw std::runtime_error("Lowest distortion eigenvalues did not converge within one percent");
        }
        if (spectrum.shear_kernel_dimension == 0)
            throw std::runtime_error("Expected nonrigid shear-null directions were lost");
        if (refinement == 3) {
            for (std::size_t mode = 0; mode < 3; ++mode) {
                double product = 0.0, left = 0.0, right = 0.0;
                for (const auto& point : section.points()) {
                    const auto a =
                        sample(*previous_section, previous_spectrum.modes[mode], point.position.x, point.position.y);
                    const auto b = sample(section, spectrum.modes[mode], point.position.x, point.position.y);
                    for (std::size_t c = 0; c < 2; ++c) {
                        product += point.weight * a[c] * b[c];
                        left += point.weight * a[c] * a[c];
                        right += point.weight * b[c] * b[c];
                    }
                }
                const double correlation = product * product / (left * right);
                std::cout << "mode=" << mode << " mesh_shape_correlation=" << correlation << '\n';
                if (correlation < 0.99)
                    throw std::runtime_error("Lowest section mode shapes did not converge");
            }
        }
        previous_section = std::make_unique<fuelsim::CrossSection>(section);
        previous_spectrum = spectrum;
        previous.clear();
        for (const auto& mode : spectrum.modes)
            previous.push_back(mode.eigenvalue);
        if (refinement == 1) {
            const auto scaled = fuelsim::test::rectangle({-0.002, 0.002},
                0.04,
                {fuelsim::test::region("scaled", 140e9, 0.3)},
                false,
                4);
            const auto scaled_spectrum = fuelsim::build_section_distortion_modes(scaled, 6);
            for (std::size_t i = 0; i < spectrum.modes.size(); ++i)
                if (std::abs(4.0 * scaled_spectrum.modes[i].eigenvalue / spectrum.modes[i].eigenvalue - 1.0) > 1.0e-8)
                    throw std::runtime_error("Distortion length scaling or elastic-modulus invariance failed");
            const auto repeated = fuelsim::build_section_distortion_modes(section, 6);
            for (std::size_t i = 0; i < repeated.modes.size(); ++i)
                if (std::abs(repeated.modes[i].eigenvalue / spectrum.modes[i].eigenvalue - 1.0) > 1.0e-12)
                    throw std::runtime_error("Repeated section eigensolve is not deterministic");
        }
    }
    // A curved isoparametric section need not have the four polynomial
    // shear-free directions present in the rectangular fixtures. Production
    // selection must use the measured kernel dimension, not assume it is four.
    std::vector<fuelsim::CartesianPoint3> curved_nodes{{-0.001, -0.01, 0.0},
        {0.001, -0.01, 0.0},
        {0.001, 0.01, 0.0},
        {-0.001, 0.01, 0.0},
        {0.0001, -0.011, 0.0},
        {0.0011, 0.001, 0.0},
        {-0.0001, 0.012, 0.0},
        {-0.0012, -0.001, 0.0}};
    fuelsim::SectionCell curved_cell{};
    curved_cell.element.nodes = {0, 1, 2, 3, 4, 5, 6, 7};
    const fuelsim::CrossSection curved(std::move(curved_nodes),
        {curved_cell},
        {fuelsim::test::region("curved", 70e9, 0.3)});
    const auto curved_modes = fuelsim::build_section_distortion_modes(curved, 3, 100);
    std::cout << "curved_section_shear_kernel=" << curved_modes.shear_kernel_dimension << '\n';
    if (curved_modes.shear_free_modes.size() != curved_modes.shear_kernel_dimension
        || curved_modes.shear_kernel_dimension >= 4)
        throw std::runtime_error("Curved section did not exercise the geometry-dependent shear kernel");
}
} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "Section distortion internal contracts");
        std::cout << std::setprecision(16);
        verify();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
