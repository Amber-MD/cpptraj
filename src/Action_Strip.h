#ifndef INC_ACTION_STRIP_H
#define INC_ACTION_STRIP_H
#include "Action.h"
// Class: Strip
/// Used to remove atoms from the state.
class Strip: public Action {
    AmberParm *oldParm;
    AmberParm *newParm;
    Frame newFrame;
    char *prefix;
    AtomMask M1;
    bool removeBoxInfo;
  public:
    Strip();
    ~Strip();

    int init();
    int setup();
    int action();
};
// Class: Unstrip
/// Signals to ActionList that the original traj parm should be restored.
class Unstrip: public Action {
  public:
    Unstrip() {}
    ~Unstrip() {}

    int setup() {return 2;}
    int action() {return 2;}
};
#endif  
