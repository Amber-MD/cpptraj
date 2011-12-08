#ifndef INC_ACTIONS_PAIRWISE_H
#define INC_ACTIONS_PAIRWISE_H
#include <vector>
#include "Action.h"
// Class: Pairwise 
/// Action to calculate nonbonded energy between pairs of atoms.
/** Functions in two ways:
  * - Calculate Lennard-Jones and Coulomb energy for each frame. 
  *   Also calculate the cumulative LJ and Coulomb energy on each atom.
  * - Calculate the Lennard-Jones and Coulomb energy between each
  *   pair of atoms in a reference structure. Calculate the difference
  *   in each pair from frame to reference (d = Ref - Frame). 
  */
class Pairwise: public Action {
    enum PairCalcType { SET_REF, COMPARE_REF, NORMAL };
    PairCalcType nb_calcType;       ///< Type of nonbonded calc being performed
    AtomMask Mask0;                 ///< Calculate energy for atoms in mask
    AtomMask RefMask;               ///< Reference mask
    AmberParm *RefParm;             ///< Reference parm
    Frame *RefFrame;                ///< Reference coordinates
    int N_ref_interactions;         ///< Number of interactions in Ref w/exclusions
    double kes;                     ///< Electrostatic constant, 1.0 when using Amber units
    DataSet *ds_vdw;                ///< Evdw dataset
    DataSet *ds_elec;               ///< Eelec dataset
    double ELJ, Eelec;              ///< Total Evdw and Eelec over all atoms
    double cut_evdw, cut_evdw1;     ///< Evdw cutoff ( Evdw < cutevdw1 && Evdw > cutevdw )
    std::vector<double> atom_evdw;  ///< Cumulative Evdw on each atom
    double cut_eelec, cut_eelec1;   ///< Eelec cutoff ( Eelec < cuteelec1 && Eelec > cuteelec )
    std::vector<double> atom_eelec; ///< Cumulative Eelec on each atom
    char *cutout;                   ///< Mol2 file prefix for atoms satisfying cutoffs
    CpptrajFile Eout;               ///< Output file for atom energies.
    /// HACK: Hold charges * 18.2223
    std::vector<double> atom_charge;
    /// Hold the exclusion list for each atom
    std::vector< std::vector<int> > exclusionList;
    /// Hold nonbond energy for a given atom pair
    struct NonbondEnergyType {
      double evdw;
      double eelec;
    };
    /// Hold nonbond energy for each pair of atoms in reference
    std::vector<NonbondEnergyType> ref_nonbondEnergy;
    /// Hold cumulative LJ and elec energy for each atom
    //std::vector<NonbondEnergyType> atom_nonbondEnergy;

    /// Set up nonbondParm for given Parm and atoms in mask
    int SetupNonbondParm(AtomMask &, AmberParm *);
    /// Calculate nonbond energy using nonbondParm for given frame
    void NonbondEnergy(Frame *, AmberParm *, AtomMask &);
    int WriteCutFrame(AmberParm *, AtomMask *, double *, Frame *, char*);
    void PrintCutAtoms(Frame *);

  public:
    Pairwise();
    ~Pairwise();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
