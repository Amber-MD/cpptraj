#ifndef INC_ACTION_DISTRMSD_H
#define INC_ACTION_DISTRMSD_H
#include "Action.h"
#include "ReferenceAction.h"
// Class: Action_DistRmsd
/// Action to calculate the distance RMSD between frame and a reference frame.
class Action_DistRmsd: public Action {
  public:
    Action_DistRmsd();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_DistRmsd(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    ReferenceAction refHolder_;
    DataSet *drmsd_;    ///< DRMSD DataSet
    AtomMask TgtMask_;  ///< Target mask.
    Frame SelectedTgt_; ///< Hold only target coords selected by TgtMask
};
#endif
