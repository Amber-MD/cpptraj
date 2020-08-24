#ifndef INC_EXEC_GRAFT_H
#define INC_EXEC_GRAFT_H
#include "Exec.h"
/// Graft part of one COORDS to another COORDS 
class Exec_Graft : public Exec {
  public:
    Exec_Graft() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Graft(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
