#ifndef INC_ENERGY_ENERGYDECOMPOSER_H
#define INC_ENERGY_ENERGYDECOMPOSER_H
#include "../CharMask.h"
namespace Cpptraj {
namespace Energy {
/// Used to break down pairwise-additive energy by atom or residue
class EnergyDecomposer {
  public:
    /// CONSTRUCTOR
    EnergyDecomposer();
  private:
    CharMask selectedAtoms_; ///< Mask of atoms that energy will be recorded for.
};
}
}
#endif
