#ifndef INC_ENERGY_ENERGYDECOMPOSER_H
#define INC_ENERGY_ENERGYDECOMPOSER_H
#include <vector>
#include "../CharMask.h"
#include "../OnlineVarT.h"
class ArgList;
class DataSet;
class DataSetList;
namespace Cpptraj {
namespace Energy {
/// Used to break down pairwise-additive energy by atom or residue
class EnergyDecomposer {
  public:
    /// CONSTRUCTOR
    EnergyDecomposer();
    /// Initialize with arguments
    int InitDecomposer(ArgList&, DataSetList&, int);
    /// Print options to stdout
    void PrintOpts() const;
  private:
    typedef std::vector< Stats<double> > EneArrayType;

    CharMask selectedAtoms_; ///< Mask of atoms that energy will be recorded for.
    DataSet* eneOut_;        ///< Will hold the average energy of each selected entity for output.
    EneArrayType energies_;  ///< Used to accumulate the average energy of each selected entity.
    int debug_;              ///< Debug level
};
}
}
#endif
