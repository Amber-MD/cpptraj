#ifndef INC_EXEC_SOURCE_H
#define INC_EXEC_SOURCE_H
#include "Exec.h"
/// Replicate LEaPs source command 
class Exec_Source : public Exec {
  public:
    Exec_Source() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Source(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
