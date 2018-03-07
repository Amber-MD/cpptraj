#ifndef INC_EXEC_RUNANALYSIS_H
#define INC_EXEC_RUNANALYSIS_H
#include "Exec.h"
class Exec_RunAnalysis : public Exec {
  public:
    Exec_RunAnalysis() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_RunAnalysis(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int DoRunAnalysis(CpptrajState&, ArgList&) const;
};
#endif
