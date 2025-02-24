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

    /// Set up an energy component output data set
    DataSet* addCompSet(DataSetList&, std::string const&);
    /// Set up selected bonds
    int setupBonds(BndArrayType const&);
    /// Set up selected angles
    int setupAngles(AngArrayType const&);
    /// Set up selected dihedrals
    int setupDihedrals(DihArrayType const&);

    /// Save energy contribution for atom if it is selected
    inline void saveEne(int, double, Darray&);
    /// Calculate bond energies
    void calcBonds(Frame const&);
    /// Calculate angle energies
    void calcAngles(Frame const&);
    /// Calculate dihedral energies
    void calcDihedrals(Frame const&);
    /// Put final average energies in output data set
    void populateOutputData(DataSet*, EneArrayType const&) const;

    CharMask selectedAtoms_; ///< Mask of atoms that energy will be recorded for.
    DataSet* eneOut_;        ///< Will hold the average energy of each selected entity for output.
    DataSet* eBndOut_;       ///< Will hold average bond energy of each selected atom for output.
    DataSet* eAngOut_;       ///< Will hold average angle energy of each selected atom for output.
    DataSet* eDihOut_;       ///< Will hold average dihedral energy of each selected atom for output.
    DataSet* eV14Out_;       ///< Will hold average 1-4 vdW energy of each selected atom for output.
    DataSet* eE14Out_;       ///< Will hold average 1-4 elec. energy of each selected atom for output.
    DataSet* eEleOut_;       ///< Will hold average elec. energy of each selected atom for output.
    DataSet* eVdwOut_;       ///< Will hold average vdW energy of each selected atom for output.
    DataFile* outfile_;      ///< Output file
    int debug_;              ///< Debug level
    bool saveComponents_;    ///< If true, save per-atom energies for individual energy components 
    EwaldOptions ewaldOpts_; ///< Hold Ewald options
    Cpptraj::Energy::Ecalc_Nonbond::CalcType nbcalctype_;

    BndArrayType bonds_;         ///< Hold all bonds to be calculated
    AngArrayType angles_;        ///< Hold all angles to be calculated
    DihArrayType dihedrals_;     ///< Hold all dihedrals to be calculated
    //AtomMask mask_;              ///< Atom mask for nonbonded calculations
    EneArrayType energies_;      ///< Used to accumulate the average energy of each selected entity.
    EneArrayType eBonds_;        ///< Accumulate average bond energy
    EneArrayType eAngles_;       ///< Accumulate average angle energy
    EneArrayType eDihedrals_;    ///< Accumulate average dihedral energy
    EneArrayType eVDW14_;        ///< Accumulate average 1-4 vdW energy
    EneArrayType eELE14_;        ///< Accumulate average 1-4 elec. energy
    EneArrayType eElec_;         ///< Accumulate average electrostatic energy
    EneArrayType eVdw_;          ///< Accumulate average vdW energy
    Topology const* currentTop_; ///< Current topology from Setup

    Ecalc_Nonbond NB_;  ///< For calculating pairwise decomposed nonbond energies

    Darray currentEne_;      ///< Hold the total energy of each atom for the current frame
    Darray currentBnd_;      ///< Hold the bond energy of each atom for the current frame
    Darray currentAng_;      ///< Hold the angle energy of each atom for the current frame
    Darray currentDih_;      ///< Hold the dihedral energy of each atom for the current frame
    Darray currentV14_;      ///< Hold the 1-4 vdW energy of each atom for the current frame
    Darray currentE14_;      ///< Hold the 1-4 elec. energy of each atom for the current frame
    Darray currentELE_;      ///< Hold elec. energy of each atom for the current frame
    Darray currentVDW_;      ///< Hold vdW energy of each atom for the current frame

    Timer t_total_; ///< Total calc time
};
}
}
#endif
