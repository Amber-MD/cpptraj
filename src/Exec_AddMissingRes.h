#ifndef INC_EXEC_ADDMISSINGRES_H
#define INC_EXEC_ADDMISSINGRES_H
#include "Exec.h"
/// <Enter description of Exec_AddMissingRes here>
class Exec_AddMissingRes : public Exec {
  public:
    Exec_AddMissingRes() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_AddMissingRes(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
