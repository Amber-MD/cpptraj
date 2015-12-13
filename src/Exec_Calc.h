#ifndef INC_EXEC_CALC_H
#define INC_EXEC_CALC_H
#include "Exec.h"
class Exec_Calc : public Exec {
  public:
    Exec_Calc() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Calc(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
