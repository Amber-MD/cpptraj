#ifndef ACTION_BOX_H
#define ACTION_BOX_H
#include "Action.h"
/// Manipulate box coords
class Action_Box : public Action {
  public:
    Action_Box();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Box(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum ModeType { SET = 0, REMOVE, AUTO };

    CoordinateInfo cInfo_; ///< For holding modified coordinate info.
    Box box_;              ///< Hold box info to be set for SET.
    ModeType mode_;        ///< How box info will be assigned.
};
#endif
