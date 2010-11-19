#ifndef INC_ACTION_DISTANCE_H
#define INC_ACTION_DISTANCE_H

// Should automatically include AmberParm.h from Action.h
#include "Action.h"
//#include "AtomMask.h" // NOTE: Should be included in Action.h instead?

class Distance: public Action {
    DataSet *dist;
    bool noimage,useMass;
    enum ImageType {NONE, ORTHO, NONORTHO};
    ImageType imageType; 
    AtomMask Mask1, Mask2;
  public:
    Distance();
    ~Distance();

    int init();
    int setup();
    int action();
};
#endif  
