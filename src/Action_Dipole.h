#ifndef INC_ACTION_DIPOLE_H
#define INC_ACTION_DIPOLE_H
class DataSet_GridFlt;
#include "Action.h"
#include "GridAction.h"
#include "CharMask.h"
class Action_Dipole : public Action, private GridAction {
  public:
    Action_Dipole();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Dipole(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print();

    DataSet_GridFlt* grid_;
    std::vector<Vec3> dipole_;
    CpptrajFile* outfile_;
    CharMask mask_;
    double max_;
    Topology* CurrentParm_;
};
#endif
