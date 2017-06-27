#ifndef INC_ACTION_UNIMAGE_H
#define INC_ACTION_UNIMAGE_H
#include "Action.h"
#include "ImagedAction.h"
/// <Enter description of Action_Unimage here>
class Action_Unimage : public Action {
  public:
    //Action_Unimage() : natoms_(0), useCenter_(false) {}
    Action_Unimage() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Unimage(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    ImagedAction image_ ; ///< Imaging routines
    CharMask mask_;
    BondArray bonds_;
    Vec3 boxCenter_; ///< Box center for current frame
    Matrix_3x3 ucell_; ///< Unit cell matrix for current frame
    Matrix_3x3 recip_; ///< Recip (frac) matrix for current frame
    Topology* CurrentParm_;
};
#endif
