#ifndef INC_ACTION_ROTATE_H
#define INC_ACTION_ROTATE_H
#include "Action.h"
class Action_Rotate : public Action {
  public:
    Action_Rotate();
  private:
    double RotMatrix_[9];
    AtomMask mask_;

    int init();
    int setup();
    int action();
    void print() {}
};
#endif
