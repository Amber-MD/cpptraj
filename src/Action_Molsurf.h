#ifndef INC_ACTION_MOLSURF_H
#define INC_ACTION_MOLSURF_H
#include "Action.h"
// Class: Molsurf
/// Wrapper for the molsurf routine in molsurf.c 
class Molsurf: public Action {
    DataSet *sasa;
    AtomMask Mask1;
  public:
    Molsurf();
    ~Molsurf();

    int init();
    int setup();
    int action();
};
#endif
