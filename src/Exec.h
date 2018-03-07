#ifndef INC_EXEC_H
#define INC_EXEC_H
#include "DispatchObject.h"
#include "CpptrajState.h"
/// An Exec is executed immediately instead of being queued.
class Exec : public DispatchObject {
  public:
    /// CONSTRUCTOR
    Exec() {}
    /// CONSTRUCTOR - Take DispatchObject::Otype
    Exec(Otype o) : DispatchObject( o ) {}
    /// Alias for CpptrajState::RetType
    typedef CpptrajState::RetType RetType; 
    /// Actually execute the command.
    virtual RetType Execute(CpptrajState&, ArgList&) = 0;
};
#endif
