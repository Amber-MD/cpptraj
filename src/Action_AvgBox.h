#ifndef INC_ACTION_AVGBOX_H
#define INC_ACTION_AVGBOX_H
#include "Action.h"
#include "OnlineVarT.h"
/// Calculate average box size 
class Action_AvgBox : public Action {
  public:
    Action_AvgBox() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AvgBox(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Stats<double> avgbox_[9]; ///< For averaging box unit cell vectors
    DataSet* boxMatrix_; ///< Hold box matrix data set
};
#endif
