#ifndef INC_ACTION_ANGLE_H
#define INC_ACTION_ANGLE_H
#include "Action.h"
// Class: Angle
/// Calculate angle between atom(s) in 3 masks
class Angle: public Action {
    DataSet *ang;
    AtomMask Mask1, Mask2, Mask3;
  public:
    Angle();

    int init();
    int setup();
    int action();
};
#endif
