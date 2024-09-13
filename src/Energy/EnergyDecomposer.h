#ifndef INC_ENERGY_ENERGYDECOMPOSER_H
#define INC_ENERGY_ENERGYDECOMPOSER_H
#include <vector>
#include "../CharMask.h"
#include "../OnlineVarT.h"
class AngleType;
class ArgList;
class BondType;
class DataFile;
class DataFileList;
class DataSet;
class DataSetList;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
/// Used to break down pairwise-additive energy by atom
class EnergyDecomposer {
  public:
    /// CONSTRUCTOR
    EnergyDecomposer();
    /// Initialize with arguments
    int InitDecomposer(ArgList&, DataSetList&, DataFileList&, int);
    /// Print options to stdout
    void PrintOpts() const;
    /// Topology-based setup
    int SetupDecomposer(Topology const&);
    /// Calculate and decompose energies for given frame.
    int CalcEne(Frame const&);
    /// Finish the calculation by putting energies in output DataSet
    int FinishCalc();
  private:
    typedef std::vector<double> Darray;
    typedef std::vector< Stats<double> > EneArrayType;
    typedef std::vector<BondType> BndArrayType;
    typedef std::vector<AngleType> AngArrayType;

    /// Set up selected bonds
    int setupBonds(BndArrayType const&);
    /// Set up selected angles
    int setupAngles(AngArrayType const&);

    /// Save energy contribution for atom if it is selected
    inline void saveEne(int, double);
    /// Calculate bond energies
    void calcBonds(Frame const&);
    /// Calculate angle energies
    void calcAngles(Frame const&);

    CharMask selectedAtoms_; ///< Mask of atoms that energy will be recorded for.
    DataSet* eneOut_;        ///< Will hold the average energy of each selected entity for output.
    DataFile* outfile_;      ///< Output file
    int debug_;              ///< Debug level

    BndArrayType bonds_;     ///< Hold all bonds to be calculated
    AngArrayType angles_;    ///< Hold all angles to be calculated
    EneArrayType energies_;  ///< Used to accumulate the average energy of each selected entity.
    Topology const* currentTop_;

    Darray currentEne_;      ///< Hold the total energy of each atom for the current frame
};
}
}
#endif
