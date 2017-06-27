#ifndef INC_ACTION_FIXIMAGEDBONDS_H
#define INC_ACTION_FIXIMAGEDBONDS_H
#include "Action.h"
#include "ImagedAction.h"
/// <Enter description of Action_FixImagedBonds here>
class Action_FixImagedBonds : public Action {
  public:
    Action_FixImagedBonds() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_FixImagedBonds(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    ImagedAction image_ ; ///< Imaging routines
    CharMask mask_;
    Vec3 boxCenter_; ///< Box center for current frame
    Matrix_3x3 ucell_; ///< Unit cell matrix for current frame
    Matrix_3x3 recip_; ///< Recip (frac) matrix for current frame
    Topology* CurrentParm_;
    std::vector<bool> atomVisited_;
    int firstSelected_;
    unsigned int Natoms_;
};
#endif
