#ifndef INC_EXEC_PARMWRITE_H
#define INC_EXEC_PARMWRITE_H
#include "Exec.h"
/// Write Topology to file.
class Exec_ParmWrite : public Exec {
  public:
    Exec_ParmWrite() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParmWrite(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
