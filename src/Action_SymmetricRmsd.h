#ifndef INC_SYMMETRICRMSD_H
#define INC_SYMMETRICRMSD_H
#include "Action.h"
#include "ReferenceAction.h"
#include "SymmetricRmsdCalc.h"
/// Action to calculate symmetry-corrected RMSD
class Action_SymmetricRmsd : public Action {
  public:
    Action_SymmetricRmsd();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_SymmetricRmsd(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}
#   ifdef MPI
    int ParallelActionInit(Parallel::Comm const& c) { return REF_.SetTrajComm( c ); }
#   endif
    DataSet* rmsd_;           ///< Output DataSet
    Action::RetType action_return_; ///< Used to indicate if frame has been modified.
    ReferenceAction REF_;     ///< Hold reference frame/traj/options
    SymmetricRmsdCalc SRMSD_; ///< Symmetric RMSD calculation.
    AtomMask tgtMask_;        ///< Atom mask selecting target atoms.
    Frame selectedTgt_;       ///< Frame holding selected target atoms.
    bool remap_;              ///< If true, re-map target frame.
    Frame remapFrame_;        ///< Original target frame re-mapped for symmetry
    SymmetricRmsdCalc::Iarray targetMap_; ///< targetMap_[oldTgt] = newTgt
};
#endif
