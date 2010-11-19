#ifndef INC_ACTION_ANGLE_H
#define INC_ACTION_ANGLE_H
// Angle
#include "Action.h"

class Angle: public Action {
    DataSet *ang;
    AtomMask Mask1, Mask2, Mask3;
  public:
    Angle();
    ~Angle();

    int init();
    int setup();
    int action();
};
#endif
