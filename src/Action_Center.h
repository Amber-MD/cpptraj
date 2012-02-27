#ifndef INC_ACTION_CENTER_H
#define INC_ACTION_CENTER_H
/// Class: Center
/// Action to center coordinates to coord origin or box origin.
#include "Action.h"
class Center: public Action {
    AtomMask Mask1;
    double box[3];
    bool origin;
  public:
    Center();

    int init();
    int setup();
    int action();
};
#endif
