#ifndef INC_ENERGY_ENERGYDECOMPOSER_H
#define INC_ENERGY_ENERGYDECOMPOSER_H
#include "../CharMask.h"
class ArgList;
namespace Cpptraj {
namespace Energy {
/// Used to break down pairwise-additive energy by atom or residue
class EnergyDecomposer {
  public:
    /// CONSTRUCTOR
    EnergyDecomposer();
    /// Initialize with arguments
    int InitDecomposer(ArgList&, int);
    /// Print options to stdout
    void PrintOpts() const;
  private:
    CharMask selectedAtoms_; ///< Mask of atoms that energy will be recorded for.
    int debug_;              ///< Debug level
};
}
}
#endif
