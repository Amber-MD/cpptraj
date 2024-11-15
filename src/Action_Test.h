#ifndef INC_ACTION_TEST_H
#define INC_ACTION_TEST_H
#include "Action.h"
#include "Timer.h"
/// <Enter description of Action_Test here>
class Action_Test : public Action {
  public:
    Action_Test() { SetHidden(true); }
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Test(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    AtomMask mask1_, mask2_;
    std::vector<int> firstAtoms_;
    //DataSet* D1_;
    //DataSet* D2_;
    CpptrajFile* outfile_;
    Timer t1_;
    Timer t2_;
    Timer t3_;
};
#endif
