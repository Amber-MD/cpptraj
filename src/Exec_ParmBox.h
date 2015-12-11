#ifndef INC_EXEC_PARMBOX_H
#define INC_EXEC_PARMBOX_H
#include "Exec.h"
/// Modify specified parm box info.
class Exec_ParmBox : public Exec {
  public:
    Exec_ParmBox() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParmBox(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
