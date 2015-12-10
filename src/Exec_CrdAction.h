#ifndef INC_EXEC_CRDACTION_H
#define INC_EXEC_CRDACTION_H
#include "Exec.h"
/// Perform action on given COORDS DataSet
class Exec_CrdAction : public Exec {
  public:
    Exec_CrdAction() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CrdAction(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    RetType DoCrdAction(CpptrajState&, ArgList&, DataSet_Coords*,
                        Action*, TrajFrameCounter const&) const;
};
#endif
