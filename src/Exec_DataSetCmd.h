#ifndef INC_EXEC_DATASETCMD_H
#define INC_EXEC_DATASETCMD_H
#include "Exec.h"
/// Process DataSet-specific command
class Exec_DataSetCmd : public Exec {
  public:
    Exec_DataSetCmd() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_DataSetCmd(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
