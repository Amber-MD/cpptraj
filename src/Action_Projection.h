#ifndef INC_ACTION_PROJECTION_H
#define INC_ACTION_PROJECTION_H
#include "Action.h"
#include "DataSet_Modes.h"
#include "ActionFrameCounter.h"
#include "Array1D.h"
/// project snapshots on normal modes
class Action_Projection : public Action, ActionFrameCounter {
  public:
    Action_Projection();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Projection(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    typedef std::vector<DataSet*> Darray;
    Darray project_;
    DataSet_Modes* modinfo_;
    int beg_;
    int end_;
    std::vector<double> sqrtmasses_;
    AtomMask mask_;
    Array1D DihedralSets_;
};
#endif
