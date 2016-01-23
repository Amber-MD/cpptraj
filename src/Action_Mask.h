#ifndef INC_ACTION_MASK_H
#define INC_ACTION_MASK_H
#include "Action.h"
#include "Trajout_Single.h"
/// Print out all atoms selected by a mask for each frame.
/** This allows use of distance-dependent masks. This does NOT modify the
  * frame or parm. 
  */
class Action_Mask: public Action {
  public:
    Action_Mask();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Mask(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    // NOTE: not setting parallel comm for output traj since only PDB/mol2 multi right now
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print() {}

    CharMask Mask1_;         ///< Atoms which will be selected each frame
    CpptrajFile* outfile_;   ///< File to write selected atom info to
    DataSet* fnum_;          ///< Hold frame numbers for selections
    DataSet* anum_;          ///< Hold selected atom numbers
    DataSet* aname_;         ///< Hold selected atom names
    DataSet* rnum_;          ///< Hold selected residue numbers
    DataSet* rname_;         ///< Hold selected residue names
    DataSet* mnum_;          ///< Hold selected molecule numbers
    int idx_;                ///< Index into data sets
    Trajout_Single outtraj_; ///< Output PDB/Mol2
    Topology* CurrentParm_;  ///< Pointer to current topology.
    CoordinateInfo cInfo_;   ///< Current coordinate info
    int debug_;
    bool writeTraj_;
};
#endif  
