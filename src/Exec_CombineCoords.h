#ifndef INC_EXEC_COMBINECOORDS_H
#define INC_EXEC_COMBINECOORDS_H
#include "Exec.h"
class Exec_CombineCoords : public Exec {
  public:
    Exec_CombineCoords() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CombineCoords(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
