#ifndef INC_ACTION_RADGYR_H
#define INC_ACTION_RADGYR_H
#include "Action.h"
/// Action to calculate the radius of gyration of atoms within a mask.
class Action_Radgyr: public Action {
  public:
    Action_Radgyr();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Radgyr(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet* rog_;
    DataSet* rogmax_;
    DataSet* rogtensor_;
    AtomMask Mask1_;
    bool calcRogmax_;
    bool calcTensor_;
    bool useMass_;
};
#endif
