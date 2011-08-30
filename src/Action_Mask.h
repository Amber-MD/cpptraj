#ifndef INC_ACTION_MASK_H
#define INC_ACTION_MASK_H
/// Class: ActionMask
/// Action that will print out all atoms selected by a mask for each frame.
/// This allows use of distance-dependent masks. This does NOT modify the
/// frame or parm. 
#include "Action.h"
class ActionMask: public Action {
    AtomMask Mask1;
    PtrajFile outfile;
    char *maskpdb;
    ArgList maskpdbarg;
  public:
    ActionMask();
    ~ActionMask();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
