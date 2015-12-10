#ifndef INC_EXEC_COMMANDS_H
#define INC_EXEC_COMMANDS_H
/*! \file Exec_Commands.h
    \brief Exec classes for very simple commands.
 */
#include "Exec.h"
/// Run queued Actions/Analyses/output Trajectories in the State
class Exec_Run : public Exec {
  public:
    Exec_Run() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Run(); }
    RetType Execute(CpptrajState& State, ArgList&) { return (CpptrajState::RetType)State.Run(); }
};
#endif
