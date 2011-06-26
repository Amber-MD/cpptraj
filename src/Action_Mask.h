#ifndef INC_ACTION_MASK_H
#define INC_ACTION_MASK_H

// Should automatically include AmberParm.h from Action.h
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
