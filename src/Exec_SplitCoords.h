#ifndef INC_EXEC_SPLITCOORDS_H
#define INC_EXEC_SPLITCOORDS_H
#include "Exec.h"
/// Split a trajectory up in various ways 
class Exec_SplitCoords : public Exec {
  public:
    Exec_SplitCoords() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SplitCoords(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
