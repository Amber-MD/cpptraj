#ifndef INC_EXEC_FLATTEN_H
#define INC_EXEC_FLATTEN_H
#include "Exec.h"
/// Flatten a 2D matrix to a 1D data set 
class Exec_Flatten : public Exec {
  public:
    Exec_Flatten() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Flatten(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
