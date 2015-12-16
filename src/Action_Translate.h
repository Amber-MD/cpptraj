#ifndef INC_ACTION_TRANSLATE_H
#define INC_ACTION_TRANSLATE_H
#include "Action.h"
class Action_Translate : public Action {
  public:
    Action_Translate();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Translate(); }
    void Help() const;
  private:
    Vec3 Trans_;
    AtomMask mask_;

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}
};
#endif
