#ifndef INC_ACTION_KEEP_H
#define INC_ACTION_KEEP_H
#include "Action.h"
class DataSet_string;
/// Keep only specified parts of the system 
class Action_Keep : public Action {
  public:
    Action_Keep();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Keep(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet_string* bridgeData_; ///< Bridging water ID data set
    int nbridge_;
};
#endif
