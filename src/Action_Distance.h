#ifndef INC_ACTION_DISTANCE_H
#define INC_ACTION_DISTANCE_H
#include "Action.h"
// Class: Distance
/// Action to calculate a distance between atoms in two masks.
class Distance: public Action {
    DataSet *dist;
    AtomMask Mask1, Mask2;
  public:
    Distance();

    int init();
    int setup();
    int action();
};
#endif  
