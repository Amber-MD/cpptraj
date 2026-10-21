#ifndef INC_EXEC_FIXATOMNAMES_H
#define INC_EXEC_FIXATOMNAMES_H
#include "Exec.h"
/// <Enter description of Exec_FixAtomNames here>
class Exec_FixAtomNames : public Exec {
  public:
    Exec_FixAtomNames() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_FixAtomNames(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
