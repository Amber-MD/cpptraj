#ifndef INC_ENERGY_ENERGYDECOMPOSER_H
#define INC_ENERGY_ENERGYDECOMPOSER_H
#include <vector>
#include "Ecalc_Nonbond.h"
#include "../CharMask.h"
#include "../EwaldOptions.h"
#include "../OnlineVarT.h"
#include "../Timer.h"
#ifdef MPI
# include "../Parallel.h"
#endif
class AngleType;
class ArgList;
class BondType;
class DataFile;
class DataFileList;
class DataSet;
class DataSetList;
class DihedralType;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
/// Used to break down pairwise-additive energy by atom
/** Comple with -DCPPTRAJ_DEBUG_ENEDECOMP for more details on the individual
  * contributions.
  */
class EnergyDecomposer {
  public:
    /// CONSTRUCTOR
    EnergyDecomposer();
    /// Print help text to stdout
    static void HelpText();
    /// Initialize with arguments
    int InitDecomposer(ArgList&, DataSetList&, DataFileList&, int);
    /// Print options to stdout
    void PrintOpts() const;
    /// Topology-based setup. Box only needed for PME
    int SetupDecomposer(Topology const&, Box const& boxIn);
    /// Calculate and decompose energies for given frame.
    int CalcEne(Frame const&);
    /// Finish the calculation by putting energies in output DataSet
    int FinishCalc();
#   ifdef MPI
    /// Reduce the decomposed array to the master rank
    int ReduceToMaster(Parallel::Comm const&);
#   endif
  private:
    typedef std::vector<double> Darray;
    typedef std::vector< Stats<double> > EneArrayType;
    typedef std::vector<BondType> BndArrayType;
    typedef std::vector<AngleType> AngArrayType;
    typedef std::vector<DihedralType> DihArrayType;

    /// Set up selected bonds
    int setupBonds(BndArrayType const&);
    /// Set up selected angles
    int setupAngles(AngArrayType const&);
    /// Set up selected dihedrals
    int setupDihedrals(DihArrayType const&);

    /// Save energy contribution for atom if it is selected
    inline void saveEne(int, double);
    /// Calculate bond energies
    void calcBonds(Frame const&);
    /// Calculate angle energies
    void calcAngles(Frame const&);
    /// Calculate dihedral energies
    void calcDihedrals(Frame const&);

    CharMask selectedAtoms_; ///< Mask of atoms that energy will be recorded for.
    DataSet* eneOut_;        ///< Will hold the average energy of each selected entity for output.
    DataFile* outfile_;      ///< Output file
    int debug_;              ///< Debug level
    EwaldOptions ewaldOpts_; ///< Hold Ewald options
    Cpptraj::Energy::Ecalc_Nonbond::CalcType nbcalctype_;

    BndArrayType bonds_;         ///< Hold all bonds to be calculated
    AngArrayType angles_;        ///< Hold all angles to be calculated
    DihArrayType dihedrals_;     ///< Hold all dihedrals to be calculated
    //AtomMask mask_;              ///< Atom mask for nonbonded calculations
    EneArrayType energies_;      ///< Used to accumulate the average energy of each selected entity.
    Topology const* currentTop_; ///< Current topology from Setup

    Ecalc_Nonbond NB_;  ///< For calculating pairwise decomposed nonbond energies

    Darray currentEne_;      ///< Hold the total energy of each atom for the current frame

    Timer t_total_; ///< Total calc time
};
}
}
#endif
