#pragma once
#include "core/cross_section.hpp"

namespace fuelsim {
// Requires an initialized PETSc session. Offline preprocessing is local to each caller.
ClassicSectionBasis build_classic_section_basis(const CrossSection& section);
} // namespace fuelsim
