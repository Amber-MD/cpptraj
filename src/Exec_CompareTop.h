#ifndef INC_EXEC_COMPARETOP_H
#define INC_EXEC_COMPARETOP_H
#include "Exec.h"
// EXPERIMENTAL ALPHA CODE
class Exec_CompareTop : public Exec {
  public:
    Exec_CompareTop() : Exec(PARM) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CompareTop(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
