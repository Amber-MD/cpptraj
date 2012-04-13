#ifndef INC_ACTION_ANGLE_H
#define INC_ACTION_ANGLE_H
#include "Action.h"
// Class: Angle
/// Calculate angle between atom(s) in 3 masks
class Angle: public Action {
  public:
    Angle();

  private:
    int init();
    int setup();
    int action();

    DataSet *ang_;
    AtomMask Mask1_;
    AtomMask Mask2_;
    AtomMask Mask3_;
};
#endif
