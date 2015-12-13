#ifndef INC_EXEC_SCALEDIHEDRALK_H
#define INC_EXEC_SCALEDIHEDRALK_H
#include "Exec.h"
/// Scale dihedral force constants in specfied parm by factor.
class Exec_ScaleDihedralK : public Exec {
  public:
    Exec_ScaleDihedralK() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ScaleDihedralK(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
