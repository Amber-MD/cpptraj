#ifndef INC_ACTION_CENTER_H
#define INC_ACTION_CENTER_H
/// Action to center coordinates to coord origin or box origin.
#include "Action.h"
class Action_Center: public Action {
  public:
    Action_Center();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Center(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum CenterMode { ORIGIN = 0, BOXCTR, REF, POINT };
    AtomMask Mask_;
    CenterMode centerMode_;
    bool useMass_;
    Vec3 refCenter_;
};
#endif
