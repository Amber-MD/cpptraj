#ifndef INC_ACTION_CONVERTTOFRAC_H
#define INC_ACTION_CONVERTTOFRAC_H
#include "Action.h"
/// <Enter description of Action_ConvertToFrac here>
class Action_ConvertToFrac : public Action {
  public:
    Action_ConvertToFrac() { SetHidden(true); }
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_ConvertToFrac(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Frame newFrame_;
};
#endif
