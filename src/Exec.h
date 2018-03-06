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
#   ifdef MPI
    /// This will determine how the Exec is executed in parallel. TODO create communicator
    enum CommEnumType {
      TRAJCOMM_MASTER = 0, ///< Only TrajComm masters execute
      ALLTHREADS           ///< All threads execute.
    };
    /// \return How exec should be executed in parallel. Default to only trajComm masters.
    virtual CommEnumType CommType() const { return TRAJCOMM_MASTER; }
#   endif
};
#endif
