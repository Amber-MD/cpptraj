#ifndef INC_EXEC_PRECISION_H
#define INC_EXEC_PRECISION_H
#include "Exec.h"
class Exec_Precision : public Exec {
  public:
    Exec_Precision() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Precision(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
