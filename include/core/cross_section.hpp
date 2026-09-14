#pragma once
#include "core/mesh.hpp"
#include "material.hpp"
#include <array>
#include <string>
#include <vector>

namespace fuelsim {
struct SectionRegion final {
    std::string name;
    IsotropicThermoelasticMaterial material;
    double temperature;
};

struct SectionCell final {
    Quad8Element element;
    std::size_t region;
};

struct SectionPoint final {
    std::array<std::size_t, 8> nodes;
    std::array<double, 8> shape;
    std::array<std::array<double, 2>, 8> gradient;
    CartesianPoint3 position;
    double weight;
    std::size_t region;
};

class CrossSection final {
  public:
    CrossSection(std::vector<CartesianPoint3> nodes,
        std::vector<SectionCell> cells,
        std::vector<SectionRegion> regions);

    const std::vector<CartesianPoint3>& nodes() const noexcept { return _nodes; }

    const std::vector<SectionPoint>& points() const noexcept { return _points; }

    const std::vector<SectionRegion>& regions() const noexcept { return _regions; }

    double area() const noexcept { return _area; }

    const CartesianPoint3& centroid() const noexcept { return _centroid; }

    const CartesianPoint3& elastic_center() const noexcept { return _elastic_center; }

  private:
    std::vector<CartesianPoint3> _nodes;
    std::vector<SectionPoint> _points;
    std::vector<SectionRegion> _regions;
    double _area = 0.0;
    CartesianPoint3 _centroid{}, _elastic_center{};
};

// Six-component arrays use tensor shear, [xx, yy, zz, xy, yz, xz].
using SectionStrain = std::array<double, 6>;

struct ClassicSectionMode final {
    enum class Kind { extension, bending_x_displacement, bending_y_displacement };
    Kind kind;
    // Field-major transverse correction per unit [axial strain, wx curvature, wy curvature].
    std::vector<double> correction;
    double auxiliary_relative_residual = 0.0;
};

struct ClassicSectionBasis final {
    CartesianPoint3 origin;
    std::array<ClassicSectionMode, 3> modes;
};

struct SectionKinematics final {
    std::array<double, 3> displacement{};
    // gradient[3*i+j] = derivative of displacement i with respect to coordinate j.
    std::array<double, 9> gradient{};
};

SectionKinematics classic_section_kinematics(const CrossSection& section,
    const ClassicSectionBasis& basis,
    std::size_t point,
    const std::array<double, 3>& q,
    const std::array<double, 3>& first,
    const std::array<double, 3>& second,
    const std::array<double, 3>& third);

struct SectionResponse final {
    std::array<double, 3> force{};
    std::array<double, 9> tangent{};
    std::vector<SectionStrain> strain, stress;
    double energy = 0.0;
};

SectionStrain classic_section_strain(const CrossSection& section,
    const ClassicSectionBasis& basis,
    std::size_t point,
    const std::array<double, 3>& generalized_strain,
    const std::array<double, 3>& generalized_strain_derivative = {});
SectionResponse evaluate_classic_section(const CrossSection& section,
    const ClassicSectionBasis& basis,
    const std::array<double, 3>& generalized_strain);
} // namespace fuelsim
