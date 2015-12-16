#ifndef INC_ACTION_ANGLE_H
#define INC_ACTION_ANGLE_H
#include "Action.h"
/// Calculate angle between atom(s) in 3 masks
class Action_Angle: public Action {
  public:
    Action_Angle();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Angle(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet *ang_;
    bool useMass_;
    AtomMask Mask1_;
    AtomMask Mask2_;
    AtomMask Mask3_;
};
#endif
