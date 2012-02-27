#ifndef INC_ACTION_MASK_H
#define INC_ACTION_MASK_H
// Class: ActionMask
/// Print out all atoms selected by a mask for each frame.
/** This allows use of distance-dependent masks. This does NOT modify the
  * frame or parm. 
  */
#include "Action.h"
class ActionMask: public Action {
    /// Atoms which will be selected each frame
    AtomMask Mask1;
    /// File to write selected atom info to
    CpptrajFile outfile;
    /// PDB output file name
    char *maskpdb;
  public:
    ActionMask();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
