#ifndef INC_ACTION_UNWRAP_H
#define INC_ACTION_UNWRAP_H
#include "Action.h"
class Action_Unwrap : public Action {
  public:
    Action_Unwrap();
  private:
    int init();
    int setup();
    int action();

    AtomMask mask_;
    Frame RefFrame_;
    Topology* RefParm_;
    bool orthogonal_;
};
#endif
