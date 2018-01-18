#ifndef INC_EXEC_UPDATEPARAMETERS_H
#define INC_EXEC_UPDATEPARAMETERS_H
#include "Exec.h"
/// <Enter description of Exec_UpdateParameters here>
class Exec_UpdateParameters : public Exec {
  public:
    Exec_UpdateParameters() : Exec(HIDDEN) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_UpdateParameters(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
