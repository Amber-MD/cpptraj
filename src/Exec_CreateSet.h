#ifndef INC_EXEC_CREATESET_H
#define INC_EXEC_CREATESET_H
#include "Exec.h"
/// <Enter description of Exec_CreateSet here>
class Exec_CreateSet : public Exec {
  public:
    Exec_CreateSet() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CreateSet(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
