#ifndef INC_EXEC_GENERATEAMBERRST_H
#define INC_EXEC_GENERATEAMBERRST_H
#include "Exec.h"
class Exec_GenerateAmberRst : public Exec {
  public:
    Exec_GenerateAmberRst() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_GenerateAmberRst(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
