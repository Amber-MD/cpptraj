#ifndef INC_EXEC_PARMBOX_H
#define INC_EXEC_PARMBOX_H
#include "Exec.h"
class Exec_ParmBox : public Exec {
  public:
    Exec_ParmBox() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParmBox(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
