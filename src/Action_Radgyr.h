#ifndef INC_ACTION_RADGYR_H
#define INC_ACTION_RADGYR_H
// Radgyr
#include "Action.h"

class Radgyr: public Action {
    DataSet *rog;
    DataSet *rogmax;
    AtomMask Mask1;
    bool useMass;
    bool calcRogmax;
  public:
    Radgyr();
    ~Radgyr();

    int init();
    int setup();
    int action();
};
#endif
