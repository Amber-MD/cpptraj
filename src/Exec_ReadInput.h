#ifndef INC_EXEC_READINPUT_H
#define INC_EXEC_READINPUT_H
#include "Exec.h"
/// Read input commands from a file.
class Exec_ReadInput : public Exec {
  public:
    Exec_ReadInput();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ReadInput(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
