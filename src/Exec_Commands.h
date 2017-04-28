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
    RetType Execute(CpptrajState& State, ArgList&) { return (RetType)State.Run(); }
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
#ifdef MPI
/// Tell CpptrajState to run parallel ensemble even with 1 thread/member
class Exec_ForceParaEnsemble : public Exec {
  public:
    Exec_ForceParaEnsemble() : Exec(HIDDEN) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ForceParaEnsemble(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
/// Exit CPPTRAJ
class Exec_Quit : public Exec {
  public:
    Exec_Quit() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Quit(); }
    RetType Execute(CpptrajState&, ArgList&) { return CpptrajState::QUIT; }
};

/// Set active reference for distance-based masks etc.
class Exec_ActiveRef : public Exec {
  public:
    Exec_ActiveRef() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ActiveRef(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.DSL().SetActiveReference( argIn );
    }
};

/// Clear data in specified lists
class Exec_Clear : public Exec {
  public:
    Exec_Clear() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Clear(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.ClearList( argIn );
    }
};

/// Remove specified data set(s)
class Exec_RemoveData : public Exec {
  public:
    Exec_RemoveData() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_RemoveData(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.RemoveDataSet( argIn );
    }
};

/// Set debug value for specified list(s)
class Exec_SetListDebug : public Exec {
  public:
    Exec_SetListDebug() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SetListDebug(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.SetListDebug( argIn );
    }
};

/// List all members of specified list(s)
class Exec_ListAll : public Exec {
  public:
    Exec_ListAll() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_ListAll(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.ListAll( argIn );
    }
};

/// Silence Actions Init/Setup output.
class Exec_SilenceActions : public Exec {
  public:
    Exec_SilenceActions() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SilenceActions(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      State.SetActionSilence( true ); return CpptrajState::OK;
    }
};

/// Process DataFile-specific command
class Exec_DataFileCmd : public Exec {
  public:
    Exec_DataFileCmd() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_DataFileCmd(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.DFL().ProcessDataFileArgs( argIn );
    }
};

/// Show results of mask expression
class Exec_SelectAtoms : public Exec {
  public:
    Exec_SelectAtoms() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SelectAtoms(); }
    RetType Execute(CpptrajState&, ArgList&);
};

/// Show results of DataSet expression
class Exec_SelectDS : public Exec {
  public:
    Exec_SelectDS() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_SelectDS(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
