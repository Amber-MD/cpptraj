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

    void print();
  private:
    int init();
    // NOTE: No setup needed for this action. Everything is done in action so 
    //       that the coords can be passed to the mask parser.
    int action();

    /// Atoms which will be selected each frame
    AtomMask Mask1_;
    /// File to write selected atom info to
    CpptrajFile outfile_;
    /// PDB output file name
    char* maskpdb_;
};
#endif  
