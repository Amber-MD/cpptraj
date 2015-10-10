#ifndef ACTION_BOX_H
#define ACTION_BOX_H
#include "Action.h"
/// Manipulate box coords
class Action_Box : public Action {
  public:
    Action_Box();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Box(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    CoordinateInfo cInfo_;
    Box box_;
    bool nobox_;
};
#endif
