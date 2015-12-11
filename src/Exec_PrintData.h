#ifndef INC_EXEC_PRINTDATA_H
#define INC_EXEC_PRINTDATA_H
#include "Exec.h"
class Exec_PrintData : public Exec {
  public:
    Exec_PrintData() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_PrintData(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
