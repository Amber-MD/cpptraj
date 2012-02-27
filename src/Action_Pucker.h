#ifndef INC_ACTION_PUCKER_H
#define INC_ACTION_PUCKER_H
// Class: Pucker
/// Calculate the ring pucker given 5 atom masks.
#include "Action.h"
class Pucker: public Action {
    DataSet *puck;
    AtomMask M1, M2, M3, M4, M5;
    int puckerMethod;
    bool amplitude;
    double offset;
    double puckermin, puckermax;
  public:
    Pucker();

    int init();
    int setup();
    int action();
};
#endif  
