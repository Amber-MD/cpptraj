#ifndef INC_ACTION_PMETEST_H
#define INC_ACTION_PMETEST_H
#include "Action.h"
#include "Energy/EwaldCalc_PME.h"
#include "Energy/EwaldCalc_LJPME.h"
#include "Ewald_ParticleMesh.h"
#include "EwaldOptions.h"
#include "Timer.h"
/// <Enter description of Action_PmeTest here>
class Action_PmeTest : public Action {
  public:
    Action_PmeTest();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_PmeTest(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    Ewald_ParticleMesh PME0_;
    Cpptraj::Energy::EwaldCalc_PME PME1_;
    Cpptraj::Energy::EwaldCalc_LJPME PME2_;
    EwaldOptions ewaldOpts_;
    AtomMask Mask1_;
    DataSet* ele_;
    DataSet* vdw_;
    Timer t_nb_;

    int method_;
    int debug_;
};
#endif