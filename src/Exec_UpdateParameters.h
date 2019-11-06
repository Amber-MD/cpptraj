#ifndef INC_EXEC_UPDATEPARAMETERS_H
#define INC_EXEC_UPDATEPARAMETERS_H
#include "Exec.h"
/// Update parameters in a topology with those from a data set. 
class Exec_UpdateParameters : public Exec {
  public:
    Exec_UpdateParameters() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_UpdateParameters(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static const char* disclaimer_;
};
#endif
