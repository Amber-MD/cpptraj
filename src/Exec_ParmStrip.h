#ifndef INC_EXEC_PARMSTRIP_H
#define INC_EXEC_PARMSTRIP_H
#include "Exec.h"
class Exec_ParmStrip : public Exec {
  public:
    Exec_ParmStrip() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParmStrip(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
