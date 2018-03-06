#ifndef INC_EXEC_CRDOUT_H
#define INC_EXEC_CRDOUT_H
#include "Exec.h"
class Exec_CrdOut : public Exec {
  public:
    Exec_CrdOut() : Exec(COORDS) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_CrdOut(); }
    RetType Execute(CpptrajState&, ArgList&);
#   ifdef MPI
    CommEnumType CommType() const { return ALLTHREADS; }
#   endif
  private:
    RetType WriteCrd(CpptrajState&, ArgList&);
#   ifdef MPI
    Parallel::Comm trajComm_;
#   endif
};
#endif
