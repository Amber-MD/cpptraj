#ifndef INC_EXEC_LOADCRD_H
#define INC_EXEC_LOADCRD_H
#include "Exec.h"
class Exec_LoadCrd : public Exec {
  public:
    Exec_LoadCrd() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_LoadCrd(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
