#ifndef INC_ACTION_DISTANCE_H
#define INC_ACTION_DISTANCE_H
#include "Action.h"
// Class: Action_Distance
/// Action to calculate a distance between atoms in two masks.
class Action_Distance: public Action {
  public:
    Action_Distance();
  private:
    DataSet *dist_;
    AtomMask Mask1_;
    AtomMask Mask2_;

    int init();
    int setup();
    int action();
};
#endif  
