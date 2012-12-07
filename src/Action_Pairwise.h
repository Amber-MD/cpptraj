#ifndef INC_ACTION_PAIRWISE_H
#define INC_ACTION_PAIRWISE_H
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
class Action_Pairwise: public Action {
  public:
    Action_Pairwise();
    ~Action_Pairwise();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Pairwise(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    enum PairCalcType { SET_REF, COMPARE_REF, NORMAL };
    PairCalcType nb_calcType_;       ///< Type of nonbonded calc being performed
    AtomMask Mask0_;                 ///< Calculate energy for atoms in mask
    AtomMask RefMask_;               ///< Reference mask
    Topology* CurrentParm_;          ///< Set to the current topology file.
    int N_ref_interactions_;         ///< Number of interactions in Ref w/exclusions
    double kes_;                     ///< Electrostatic constant, 1.0 when using Amber units
    DataSet* ds_vdw_;                ///< Evdw dataset
    DataSet* ds_elec_;               ///< Eelec dataset
    double ELJ_;                     ///< Total VDW energy over all selected atoms.
    double Eelec_;                   ///< Total elec. energy over all selected atoms.
    double cut_evdw_;                ///< Min Evdw cutoff
    double cut_evdw1_;               ///< Max Evdw cutoff
    std::vector<double> atom_evdw_;  ///< Cumulative Evdw on each atom
    double cut_eelec_;               ///< Min Eelec cutoff
    double cut_eelec1_;              ///< Max Eelec cutoff
    std::vector<double> atom_eelec_; ///< Cumulative Eelec on each atom
    std::string cutout_;             ///< Mol2 file prefix for atoms satisfying cutoffs
    CpptrajFile Eout_;               ///< Output file for atom energies.
    /// HACK: Hold charges * 18.2223
    std::vector<double> atom_charge_;
    /// Hold nonbond energy for a given atom pair
    struct NonbondEnergyType {
      double evdw;
      double eelec;
    };
    /// Hold nonbond energy for each pair of atoms in reference
    std::vector<NonbondEnergyType> ref_nonbondEnergy_;
    /// Hold cumulative LJ and elec energy for each atom
    //std::vector<NonbondEnergyType> atom_nonbondEnergy;

    /// Set up nonbondParm for given Parm and atoms in mask
    int SetupNonbondParm(AtomMask &, Topology *);
    /// Calculate nonbond energy using nonbondParm for given frame
    void NonbondEnergy(Frame *, Topology *, AtomMask &);
    int WriteCutFrame(int, Topology const&, AtomMask const&, std::vector<double> const&, 
                      Frame const&, std::string const&);
    void PrintCutAtoms(Frame *,int);
};
#endif  
