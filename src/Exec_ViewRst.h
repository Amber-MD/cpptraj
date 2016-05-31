#ifndef INC_EXEC_VIEWRST_H
#define INC_EXEC_VIEWRST_H
#include "Exec.h"
// EXPERIMENTAL ALPHA CODE
/// View Amber restraints
class Exec_ViewRst : public Exec {
  public:
    Exec_ViewRst() : Exec(HIDDEN) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ViewRst(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
