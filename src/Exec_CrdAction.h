#ifndef INC_EXEC_CRDACTION_H
#define INC_EXEC_CRDACTION_H
#include "Exec.h"
class Exec_CrdAction : public Exec {
  public:
    Exec_CrdAction() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CrdAction(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
