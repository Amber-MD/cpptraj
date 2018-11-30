#ifndef INC_EXEC_TRAJ_H
#define INC_EXEC_TRAJ_H
/*! \file Exec_Traj.h
    \brief Exec classes for trajectory commands.
 */ 
#include "Exec.h"
/// Add input trajectory to State.
class Exec_Trajin : public Exec {
  public:
    Exec_Trajin() : Exec(TRAJ) {}
    void Help(ArgList&) const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Trajin(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.AddInputTrajectory( argIn );
    }
};
/// Add input ensemble to State.
class Exec_Ensemble : public Exec {
  public:
    Exec_Ensemble() : Exec(TRAJ) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Ensemble(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.AddInputEnsemble( argIn );
    }
};
/// Add reference coordinates to State.
class Exec_Reference : public Exec {
  public:
    Exec_Reference() : Exec(TRAJ) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Reference(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.AddReference( argIn.GetStringNext(), argIn );
    }
};

class Exec_Trajout : public Exec {
  public:
    Exec_Trajout() : Exec(TRAJ) {}
    void Help(ArgList&) const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Trajout(); }
    RetType Execute(CpptrajState& State, ArgList& argIn) {
      return (RetType)State.AddOutputTrajectory( argIn );
    }
};

/// Set expected ensemble size to improve ensemble setup efficiency in parallel.
class Exec_EnsembleSize : public Exec {
  public:
    Exec_EnsembleSize() : Exec(TRAJ) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_EnsembleSize(); }
    RetType Execute(CpptrajState&, ArgList&);
};
#endif
