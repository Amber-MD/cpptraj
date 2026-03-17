#ifndef INC_EXEC_FLUSH_H
#define INC_EXEC_FLUSH_H
#include "Exec.h"
/// <Enter description of Exec_Flush here>
class Exec_Flush : public Exec {
  public:
    Exec_Flush() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Flush(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
