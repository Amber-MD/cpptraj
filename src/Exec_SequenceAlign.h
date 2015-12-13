#ifndef INC_EXEC_SEQUENCEALIGN_H
#define INC_EXEC_SEQUENCEALIGN_H
#include "Exec.h"
// EXPERIMENTAL ALPHA CODE
class Exec_SequenceAlign : public Exec {
  public:
    Exec_SequenceAlign() : Exec(HIDDEN) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SequenceAlign(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
