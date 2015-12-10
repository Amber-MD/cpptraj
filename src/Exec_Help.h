#ifndef INC_EXEC_HELP_H
#define INC_EXEC_HELP_H
#include "Exec.h"
class Exec_Help : public Exec {
  public:
    Exec_Help() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Help(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
