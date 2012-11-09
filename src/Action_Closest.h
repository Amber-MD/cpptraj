#ifndef INC_ACTION_CLOSEST_H
#define INC_ACTION_CLOSEST_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_Closest
/// Modify the state so that only the closest solvent molecules are kept.
class Action_Closest: public Action, ImagedAction {
  public:
    Action_Closest();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Closest(); }
    static void Help();

    ~Action_Closest();

    void Print() {}
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    DataFile *outFile_;
    DataSet *framedata_;
    DataSet *moldata_;
    DataSet *distdata_;
    DataSet *atomdata_;

    int Nclosest_;          ///< Index into Closest datasets.
    std::string prefix_;          ///< Output topology prefix.
    int closestWaters_;     ///< # of waters to keep.
    bool firstAtom_;        ///< If true just calc based on solvent first atom.
    AtomMask stripMask_;    ///< Mask including solvent and kept waters.
    AtomMask distanceMask_; ///< Mask of atoms to calculate distance from solvent to.
    Topology *newParm_;     ///< New topology with solute and kept waters.
    int NsolventMolecules_; ///< # of solvent molecules in original topology.
    int debug_;
    Frame newFrame_;        ///< New frame with solute and kept waters.
    /// Hold atom #s of kept solvent in new frame.
    std::vector<int> keptWaterAtomNum_;

    /** The moldist structure is used in order to preserve the original
      * solvent molecule numbers after sorting. */
    struct MolDist {
      int mol;        ///< Original solvent molecule number (starts from 1).
      double D;       ///< Closest distance of solvent molecule to atoms in distanceMask.
      AtomMask mask;  ///< Original topology solvent molecule atom mask.
    };
    /// Return true if the first molecule is closer than the second
    struct moldist_cmp {
      inline bool operator()(MolDist first, MolDist second) const {
        return (first.D < second.D);
      }
    };
    std::vector<MolDist> SolventMols_;
};
#endif  
