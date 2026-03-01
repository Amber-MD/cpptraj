#ifndef INC_EXEC_DESC_H
#define INC_EXEC_DESC_H
#include "Exec.h"
/// <Enter description of Exec_Desc here>
class Exec_Desc : public Exec {
  public:
    Exec_Desc();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Desc(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    static int desc_atom(Topology const&, int);
};
#endif
