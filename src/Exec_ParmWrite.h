#ifndef INC_EXEC_PARMWRITE_H
#define INC_EXEC_PARMWRITE_H
#include "Exec.h"
class Exec_ParmWrite : public Exec {
  public:
    Exec_ParmWrite() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParmWrite(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
