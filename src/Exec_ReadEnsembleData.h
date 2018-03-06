#ifndef INC_EXEC_READENSEMBLEDATA_H
#define INC_EXEC_READENSEMBLEDATA_H
#include "Exec.h"
/// <Enter description of Exec_ReadEnsembleData here>
class Exec_ReadEnsembleData : public Exec {
  public:
    Exec_ReadEnsembleData() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ReadEnsembleData(); }
    RetType Execute(CpptrajState&, ArgList&);
#   ifdef MPI
    CommEnumType CommType() const { return ALLTHREADS; }
#   endif
};
#endif
