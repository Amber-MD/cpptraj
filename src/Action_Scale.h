#ifndef INC_ACTION_SCALE_H
#define INC_ACTION_SCALE_H
#include "Action.h"
class Action_Scale : public Action {
  public:
    Action_Scale();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Scale(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    AtomMask mask_;
    double sx_;
    double sy_;
    double sz_;
};
#endif
