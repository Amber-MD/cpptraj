#ifndef INC_ACTION_ANGLE_H
#define INC_ACTION_ANGLE_H
#include "Action.h"
// Class: Action_Angle
/// Calculate angle between atom(s) in 3 masks
class Action_Angle: public Action {
  public:
    Action_Angle();
    void print() {}
  private:
    int init();
    int setup();
    int action();

    DataSet *ang_;
    bool useMass_;
    AtomMask Mask1_;
    AtomMask Mask2_;
    AtomMask Mask3_;
};
#endif
