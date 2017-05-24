#ifndef INC_EXEC_CHANGE_H
#define INC_EXEC_CHANGE_H
#include "Exec.h"
/// Change things in a topology 
class Exec_Change : public Exec {
  public:
    Exec_Change() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Change(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    int ChangeResidueName(Topology&, ArgList&) const;
    int ChangeAtomName(Topology&, ArgList&) const;
};
#endif
