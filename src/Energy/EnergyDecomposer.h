#ifndef INC_ENERGY_ENERGYDECOMPOSER_H
#define INC_ENERGY_ENERGYDECOMPOSER_H
#include <vector>
#include "EwaldCalc_Decomp_PME.h"
#include "../AtomMask.h"
#include "../CharMask.h"
#include "../EwaldOptions.h"
#include "../ExclusionArray.h"
#include "../OnlineVarT.h"
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
class EnergyDecomposer {
  public:
    /// CONSTRUCTOR
    EnergyDecomposer();
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
  private:
    static const double QFAC_; ///< Coulomb prefactor

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
    /// Calculate simple nonbonded energies, no cutoff
    void calcNB_simple(Frame const&);

    CharMask selectedAtoms_; ///< Mask of atoms that energy will be recorded for.
    DataSet* eneOut_;        ///< Will hold the average energy of each selected entity for output.
    DataFile* outfile_;      ///< Output file
    int debug_;              ///< Debug level
    bool use_pme_;           ///< If true use PME for the nonbonds
    EwaldOptions ewaldOpts_; ///< Hold Ewald options

    BndArrayType bonds_;         ///< Hold all bonds to be calculated
    AngArrayType angles_;        ///< Hold all angles to be calculated
    DihArrayType dihedrals_;     ///< Hold all dihedrals to be calculated
    //AtomMask mask_;              ///< Atom mask for nonbonded calculations
    ExclusionArray Excluded_;    ///< Hold excluded atoms lists for each selected atom
    EneArrayType energies_;      ///< Used to accumulate the average energy of each selected entity.
    Topology const* currentTop_; ///< Current topology from Setup

    EwaldCalc_Decomp_PME PME_;  ///< For calculating pairwise decomposed PME energies

    Darray currentEne_;      ///< Hold the total energy of each atom for the current frame
};
}
}
#endif
