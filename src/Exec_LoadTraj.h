#ifndef INC_EXEC_LOADTRAJ_H
#define INC_EXEC_LOADTRAJ_H
#include "Exec.h"
class Exec_LoadTraj : public Exec {
  public:
    Exec_LoadTraj() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_LoadTraj(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
