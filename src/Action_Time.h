#ifndef INC_ACTION_TIME_H
#define INC_ACTION_TIME_H
#include "Action.h"
/// Add/remove/modify time information in frame. 
class Action_Time : public Action {
  public:
    Action_Time();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Time(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    double time0_;         ///< Initial time
    double dt_;            ///< Time step.
    enum ModeType { ADD = 0, MODIFY, REMOVE };
    ModeType mode_;        ///< Modification mode.
    CoordinateInfo cInfo_; ///< The modified CoordinateInfo
};
#endif
