#ifndef INC_ACTION_CLOSEST_H
#define INC_ACTION_CLOSEST_H
#include "Action.h"
#include "ImagedAction.h"
/// Modify the state so that only the closest solvent molecules are kept.
class Action_Closest: public Action {
  public:
    Action_Closest();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Closest(); }
    void Help() const;
    ~Action_Closest();
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print() {}
#   ifdef CUDA
    double* GPU_MEM_;       ///< Memory block to be sent to GPU.
    double* V_atom_coords_; ///< Hold coordinates for selected solvent atoms.
    double* U_atom_coords_; ///< Hold coordinates for selected solute atoms.
    double* V_distances_;   ///< Hold closest distance for each solvent molecule.
#   else
    typedef std::vector<double> Darray;
    Darray U_cell0_coords_; ///< Hold selected solute atoms, wrapped to cell0 if non-ortho.
#   endif

    ImagedAction image_;    ///< Imaging routines.
    DataFile *outFile_;     ///< Output file for data on closest molecules
    DataSet *framedata_;    ///< Frame number for each closest molecule.
    DataSet *moldata_;      ///< Mol# for each closest molecule.
    DataSet *distdata_;     ///< Closest distance of each molecule.
    DataSet *atomdata_;     ///< First atom of each closest molecule.
    int Nclosest_;          ///< Index into Closest molecule DataSets.
    std::string prefix_;    ///< Output topology prefix.
    std::string parmoutName_; ///< Output topology file name.
    int closestWaters_;     ///< Closest # of molecules to keep.
    int targetNclosest_;    ///< Original target # of closest molecules to keep.
    bool firstAtom_;        ///< If true just calc based on molecule first atom.
    bool useMaskCenter_;    ///< If true use geometric center of mask.
    AtomMask stripMask_;    ///< Mask including all solute and closest molecules.
    AtomMask distanceMask_; ///< Mask of atoms to calculate distance from solvent to.
    CharMask solventMask_;  ///< Optional mask selecting solvent.
    Topology *newParm_;     ///< New topology with solute and closest molecules.
    int NsolventMolecules_; ///< # of solvent molecules in SolventMols.
    int debug_;
    Frame newFrame_;        ///< New frame with solute and kept waters.
    typedef std::vector<int> Iarray;
    /// Hold atom #s of kept solvent in new frame for placing into stripMask.
    Iarray keptWaterAtomNum_;
    /** The moldist structure is used in order to preserve the original
      * solvent molecule numbers after sorting. */
    struct MolDist {
      int mol;             ///< Original solvent molecule number (starts from 1).
      double D;            ///< Closest distance of solvent molecule to atoms in distanceMask.
      AtomMask mask;       ///< Entire solvent molecule atom mask.
      Iarray solventAtoms; ///< Actual solvent atom #s to calc distance to.
    };
    /// Return true if the first molecule is closer than the second
    struct moldist_cmp {
      inline bool operator()(MolDist const& first, MolDist const& second) const {
        return (first.D < second.D);
      }
    };
    std::vector<MolDist> SolventMols_;
};
#endif  
