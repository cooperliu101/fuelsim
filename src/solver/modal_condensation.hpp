#pragma once
#include <memory>
#include <petscksp.h>
#include <vector>

namespace fuelsim {
// Private solver implementation: no PETSc types enter an installed header.
class ModalStaticCondensation final {
  public:
    ModalStaticCondensation(Mat matrix, const std::vector<bool>& eliminated);
    ~ModalStaticCondensation();
    void solve(Vec rhs, Vec solution) const;

  private:
    class Impl;
    std::unique_ptr<Impl> _impl;
};
} // namespace fuelsim
