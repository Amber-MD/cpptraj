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

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
  private:
    enum PairCalcType { SET_REF, COMPARE_REF, NORMAL };
    PairCalcType nb_calcType;       ///< Type of nonbonded calc being performed
    AtomMask Mask0;                 ///< Calculate energy for atoms in mask
    AtomMask RefMask;               ///< Reference mask
    Topology *RefParm;             ///< Reference parm
    Topology* CurrentParm_;
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
    std::string cutout_;                   ///< Mol2 file prefix for atoms satisfying cutoffs
    CpptrajFile Eout;               ///< Output file for atom energies.
    /// HACK: Hold charges * 18.2223
    std::vector<double> atom_charge;
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
    int SetupNonbondParm(AtomMask &, Topology *);
    /// Calculate nonbond energy using nonbondParm for given frame
    void NonbondEnergy(Frame *, Topology *, AtomMask &);
    int WriteCutFrame(int, Topology *, AtomMask&, std::vector<double> const&, 
                      Frame *, std::string const&);
    void PrintCutAtoms(Frame *,int);

};
#endif  
