#ifndef INC_ACTION_DISTANCE_H
#define INC_ACTION_DISTANCE_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_Distance
/// Action to calculate a distance between atoms in two masks.
class Action_Distance: public Action, ImagedAction {
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
