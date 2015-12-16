#ifndef ACTION_CREATECRD_H
#define ACTION_CREATECRD_H
#include "Action.h"
#include "DataSet_Coords_CRD.h"
class Action_CreateCrd : public Action {
  public:
    Action_CreateCrd();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_CreateCrd(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet_Coords_CRD* coords_;
    int pindex_;
    bool check_;
};
#endif
