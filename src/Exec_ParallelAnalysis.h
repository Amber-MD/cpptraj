#ifndef INC_EXEC_PARALLELANALYSIS_H
#define INC_EXEC_PARALLELANALYSIS_H
#include "Exec.h"
/// <Enter description of Exec_ParallelAnalysis here>
class Exec_ParallelAnalysis : public Exec {
  public:
    Exec_ParallelAnalysis() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ParallelAnalysis(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
