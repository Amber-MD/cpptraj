#ifndef INC_EXEC_PARMSOLVENT_H
#define INC_EXEC_PARMSOLVENT_H
#include "Exec.h"
class Exec_ParmSolvent : public Exec {
  public:
    Exec_ParmSolvent() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParmSolvent(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
