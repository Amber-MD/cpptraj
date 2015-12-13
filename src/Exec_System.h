#ifndef INC_EXEC_SYSTEM_H
#define INC_EXEC_SYSTEM_H
#include "Exec.h"
/// Run a system command
class Exec_System : public Exec {
  public:
    Exec_System() : Exec(SYSTEM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_System(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
