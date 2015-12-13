#ifndef INC_ACTION_PRINCIPAL_H
#define INC_ACTION_PRINCIPAL_H
#include "Action.h"
#include "DataSet_Mat3x3.h"
#include "DataSet_Vector.h"
class Action_Principal : public Action {
  public:
    Action_Principal();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Principal(); }
    void Help() const;
  private:
    bool doRotation_;
    bool useMass_;
    int debug_;
    AtomMask mask_;
    CpptrajFile* outfile_;
    DataSet_Mat3x3* vecData_;
    DataSet_Vector* valData_;

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}
};
#endif
