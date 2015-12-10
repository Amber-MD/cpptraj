#ifndef INC_ACTION_ROTATE_H
#define INC_ACTION_ROTATE_H
#include "Action.h"
#include "DataSet_Mat3x3.h"
class Action_Rotate : public Action {
  public:
    Action_Rotate();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Rotate(); }
    void Help() const;
  private:
    Matrix_3x3 RotMatrix_;
    AtomMask mask_;
    DataSet_Mat3x3* rmatrices_;
    bool inverse_;

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}
};
#endif
