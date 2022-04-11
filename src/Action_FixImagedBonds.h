#ifndef INC_ACTION_FIXIMAGEDBONDS_H
#define INC_ACTION_FIXIMAGEDBONDS_H
#include "Action.h"
#include "ImageOption.h"
#include "CharMask.h"
/// Fix bonds which have been distorted by e.g. by atom imaging. 
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

    ImageOption imageOpt_; ///< Used to determine if imaging should be used
    CharMask mask_;
    Vec3 boxCenter_; ///< Box center for current frame
    Topology* CurrentParm_;
    std::vector<bool> atomVisited_;
    int firstSelected_;
    unsigned int Natoms_;
};
#endif
