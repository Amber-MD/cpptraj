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

/// Tell CpptrajState to ignore errors if possible
class Exec_NoExitOnError : public Exec {
  public:
    Exec_NoExitOnError() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_NoExitOnError(); }
    RetType Execute(CpptrajState&, ArgList&);
};

/// Tell CpptrajState not to use a progress bar during Run.
class Exec_NoProgress : public Exec {
  public:
    Exec_NoProgress() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_NoProgress(); }
    RetType Execute(CpptrajState&, ArgList&);
};

/// Exit CPPTRAJ
class Exec_Quit : public Exec {
  public:
    Exec_Quit() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Quit(); }
    RetType Execute(CpptrajState&, ArgList&) { return CpptrajState::QUIT; }
};

#endif
