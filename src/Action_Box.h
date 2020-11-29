#ifndef ACTION_BOX_H
#define ACTION_BOX_H
#include "Action.h"
#include "BoxArgs.h"
/// Manipulate box coords
class Action_Box : public Action {
  public:
    Action_Box();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Box(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum ModeType { SET = 0, REMOVE, AUTO };
    enum RadiiType { UNSPECIFIED = 0, GB, PARSE, VDW, NONE };

    CoordinateInfo cInfo_; ///< For holding modified coordinate info.
    BoxArgs boxArgs_;      ///< Hold arguments for setting box (SET).
    ModeType mode_;        ///< How box info will be assigned.
    double offset_;        ///< Offset for AUTO
    RadiiType radiiMode_;  ///< Radii type to use for AUTO
    std::vector<double> Radii_; ///< Hold radius for each atom for AUTO
};
#endif
