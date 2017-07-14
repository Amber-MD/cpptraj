#ifndef INC_ACTION_LIPIDORDER_H
#define INC_ACTION_LIPIDORDER_H
#include "Action.h"
/// <Enter description of Action_LipidOrder here>
class Action_LipidOrder : public Action {
  public:
    Action_LipidOrder() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_LipidOrder(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    CharMask mask_;
};
#endif
