#ifndef INC_EXEC_DATAFILE_H
#define INC_EXEC_DATAFILE_H
#include "Exec.h"
/// Add a new DataFile to DFL with specified DataSets, to be written later.
class Exec_CreateDataFile : public Exec {
  public:
    Exec_CreateDataFile() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CreateDataFile(); }
    RetType Execute(CpptrajState&, ArgList&);
};

/// Write DataFile with specified DataSets immediately, or force write of all DataFiles in State
class Exec_WriteDataFile : public Exec {
  public:
    Exec_WriteDataFile() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_WriteDataFile(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
