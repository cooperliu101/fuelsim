#pragma once
#include "core/modal_beam.hpp"

namespace fuelsim {
ReducedSectionBasis
build_reduced_section_basis(const CrossSection& section, bool torsion, std::size_t enrichment_modes);
} // namespace fuelsim
