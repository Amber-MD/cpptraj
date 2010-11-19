#ifndef INC_ACTION_CENTER_H
#define INC_ACTION_CENTER_H
// Center
#include "Action.h"

class Center: public Action {
    AtomMask Mask1;
//    double *mass;
    double box[3];
    bool origin;
    bool useMass;
  public:
    Center();
    ~Center();

    int init();
    int setup();
    int action();
};
#endif
