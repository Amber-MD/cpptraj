#ifndef INC_ACTION_ENEDECOMP_H
#define INC_ACTION_ENEDECOMP_H
#include "Action.h"
#include "Energy/EnergyDecomposer.h"
/// Used to decompose the pairwise additive energy of a system 
class Action_EneDecomp : public Action {
  public:
    /// CONSTRUCTOR
    Action_EneDecomp() {}
    /// Allocator
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_EneDecomp(); }
    /// Help text
    void Help() const;
  private:
    /// Initialization
    Action::RetType Init(ArgList&, ActionInit&, int);
    /// Topology-based setup
    Action::RetType Setup(ActionSetup&);
    /// Do the action
    Action::RetType DoAction(int, ActionFrame&);
    /// Print results/finish calculations
    void Print();
#   ifdef MPI
    int SyncAction();
#   endif

    Cpptraj::Energy::EnergyDecomposer eneDecomp_; ///< Do the actual decomposition
#   ifdef MPI
    Parallel::Comm trajComm_; ///< Across-trajectory communicator
#   endif
};
#endif
