#ifndef INC_EXEC_READDATA_H
#define INC_EXEC_READDATA_H
#include "Exec.h"
class Exec_ReadData : public Exec {
  public:
    Exec_ReadData() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ReadData(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
