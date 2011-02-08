#ifndef INC_ACTION_PUCKER_H
#define INC_ACTION_PUCKER_H

#include "Action.h"

class Pucker: public Action {
    DataSet *puck;
    AtomMask M1, M2, M3, M4, M5;
    int puckerMethod;
    bool amplitude;
    bool useMass;
    double offset;
  public:
    Pucker();
    ~Pucker();

    int init();
    int setup();
    int action();
};
#endif  
