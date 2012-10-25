#ifndef INC_ACTION_MASK_H
#define INC_ACTION_MASK_H
#include "Action.h"
// Class: Action_Mask
/// Print out all atoms selected by a mask for each frame.
/** This allows use of distance-dependent masks. This does NOT modify the
  * frame or parm. 
  */
class Action_Mask: public Action {
  public:
    Action_Mask();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Mask(); }
    static void Help();

    void Print();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    /// Atoms which will be selected each frame
    AtomMask Mask1_;
    /// File to write selected atom info to
    CpptrajFile outfile_;
    /// PDB output file name
    std::string maskpdb_;
    Topology* CurrentParm_;
    int debug_;
};
#endif  
