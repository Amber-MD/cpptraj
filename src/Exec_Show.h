#ifndef INC_EXEC_SHOW_H
#define INC_EXEC_SHOW_H
#include "Exec.h"
/// Used to print current script variables to stdout
class Exec_Show : public Exec {
  public:
    Exec_Show() : Exec(CONTROL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Show(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
