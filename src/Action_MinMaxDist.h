#ifndef INC_ACTION_MINMAXDIST_H
#define INC_ACTION_MINMAXDIST_H
#include "Action.h"
/// Record the min/max distance between atoms/residues/molecules 
class Action_MinMaxDist : public Action {
  public:
    Action_MinMaxDist() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_MinMaxDist(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    AtomMask mask1_;
    AtomMask mask2_;
};
#endif
