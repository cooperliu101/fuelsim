#pragma once
#include "core/modal_beam.hpp"

namespace fuelsim {
ReducedSectionBasis
build_reduced_section_basis(const CrossSection& section, bool torsion, std::size_t enrichment_modes);
// Timoshenko coordinates: independent axial rotations replace -w'. Their
// Poisson corrections follow rotation gradients. Both rotations must exist.
ReducedSectionBasis make_mixed_bending_basis(ReducedSectionBasis basis);
// Use independent amplitudes for derivative fields completely represented by
// displacement modes. Retain every derivative field not contained in that span.
ReducedSectionBasis make_independent_section_basis(const CrossSection& section, ReducedSectionBasis basis);
// Conforming piecewise-linear width weights at mesh-defined y coordinates.
// Complete the seed displacement space, then retain its localized products.
ReducedSectionBasis
enrich_section_width(const CrossSection& section, ReducedSectionBasis basis, const std::vector<double>& lines);
} // namespace fuelsim
