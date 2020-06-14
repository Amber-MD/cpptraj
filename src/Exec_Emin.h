#ifndef INC_EXEC_EMIN_H
#define INC_EXEC_EMIN_H
#include "Exec.h"
/// <Enter description of Exec_Emin here>
class Exec_Emin : public Exec {
  public:
    Exec_Emin() : Exec(COORDS) { SetHidden(true); }
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Emin(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
