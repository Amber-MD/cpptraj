#ifndef INC_ACTION_BOUNDS_H
#define INC_ACTION_BOUNDS_H
#include "Action.h"
/// Report the min/max XYZ values for atoms in mask.
class Action_Bounds : public Action {
  public:
    Action_Bounds();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Bounds(); }
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
    AtomMask mask_;
    CpptrajFile* outfile_;
    double max_[3];
    double min_[3];
    Vec3 dxyz_;
    int offset_;
    DataSet* grid_;
    DataSet* ds_xmin_;
    DataSet* ds_ymin_;
    DataSet* ds_zmin_;
    DataSet* ds_xmax_;
    DataSet* ds_ymax_;
    DataSet* ds_zmax_;
};
#endif
