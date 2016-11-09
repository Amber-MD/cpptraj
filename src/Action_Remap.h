#ifndef INC_ACTION_REMAP_H
#define INC_ACTION_REMAP_H
#include "Action.h"
class Action_Remap : public Action {
  public:
    Action_Remap();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Remap(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

  std::vector<int> Map_; ///< Reference[refidx] = tgtidx
};
#endif
