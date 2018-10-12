#ifndef INC_ACTION_ALIGN_H
#define INC_ACTION_ALIGN_H
#include "Action.h"
#include "ReferenceAction.h"
/// Action to align structure onto a reference.
class Action_Align: public Action {
  public:
    Action_Align();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Align(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    ReferenceAction REF_; ///< Hold reference frame/traj/options
    AtomMask tgtMask_;    ///< Mask of selected target atoms.
    AtomMask movMask_;    ///< Mask of atoms that will be moved.
    int debug_;
    bool useMass_;        ///< If true, mass-weight the fit.
    bool moveSpecified_;  ///< If true, move mask was specified.
    Vec3 tgtTrans_;       ///< Hold translation to origin.
    Matrix_3x3 rot_;      ///< Hold best-fit rotation matrix.
    Frame tgtFrame_;      ///< Hold selected target atoms.
};
#endif
