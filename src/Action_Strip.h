#ifndef INC_ACTION_STRIP_H
#define INC_ACTION_STRIP_H

// Should automatically include AmberParm.h from Action.h
#include "Action.h"

class Strip: public Action {
    AmberParm *oldParm;
    AmberParm *newParm;
    Frame newFrame;
    char *prefix;
    AtomMask M1;
  public:
    Strip();
    ~Strip();

    int init();
    int setup();
    int action();
};

class Unstrip: public Action {
  public:
    Unstrip() {}
    ~Unstrip() {}

    int setup() {return 2;}
    int action() {return 2;}
};
#endif  
