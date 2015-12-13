#ifndef INC_ACTION_UNSTRIP_H
#define INC_ACTION_UNSTRIP_H
#include "Action.h"
/// Signals to ActionList that the original traj top/frame should be restored.
class Action_Unstrip: public Action {
  public:
    Action_Unstrip() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Unstrip(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int) { return Action::OK;                 }
    Action::RetType Setup(ActionSetup&)              { return Action::USE_ORIGINAL_FRAME; }
    Action::RetType DoAction(int, ActionFrame&)      { return Action::USE_ORIGINAL_FRAME; }
    void Print() {}
};
#endif
