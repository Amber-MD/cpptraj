#ifndef INC_ACTION_MASK_H
#define INC_ACTION_MASK_H
#include "Action.h"
#include "TrajectoryFile.h"
// Class: Action_Mask
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
    void Print() {}

    int ensembleNum_;
    CharMask Mask1_;         ///< Atoms which will be selected each frame
    CpptrajFile* outfile_;   ///< File to write selected atom info to
    std::string maskpdb_;    ///< Traj output file name
    Topology* CurrentParm_;
    CoordinateInfo currentCoordInfo_;
    int debug_;
    TrajectoryFile::TrajFormatType trajFmt_; ///< Output trajectory format
    const char* trajOpt_;    ///< Output trajectory options
};
#endif  
