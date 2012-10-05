#ifndef INC_ACTION_TRANSLATE_H
#define INC_ACTION_TRANSLATE_H
#include "Action.h"
class Action_Translate : public Action {
  public:
    Action_Translate();
  private:
    double Trans_[3];
    AtomMask mask_;

    int init();
    int setup();
    int action();
    void print() {}
};
#endif
