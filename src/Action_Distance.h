#ifndef INC_ACTION_DISTANCE_H
#define INC_ACTION_DISTANCE_H
#include "Action.h"
// Class: Distance
/// Action to calculate a distance between atoms in two masks.
class Distance: public Action {
    DataSet *dist;
    bool noimage;
    int imageType; 
    AtomMask Mask1, Mask2;
  public:
    Distance();
    ~Distance();

    int init();
    int setup();
    int action();
};
#endif  
