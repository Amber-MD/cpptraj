#ifndef INC_ACTION_CENTER_H
#define INC_ACTION_CENTER_H
/// Class: Center
/// Action to center coordinates to coord origin or box origin.
#include "Action.h"
class Center: public Action {
  public:
    Center();
  private:
    int init();
    int setup();
    int action();

    AtomMask Mask_;
    bool origin_;
};
#endif
