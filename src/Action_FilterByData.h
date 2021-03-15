#ifndef INC_ACTION_FILTERBYDATA_H
#define INC_ACTION_FILTERBYDATA_H
#include "Action.h"
#include "DataFilter.h"
/// Filter out frames by DataSet
class Action_FilterByData : public Action {
  public:
    Action_FilterByData();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_FilterByData(); }
    void Help() const;
  private:
    // Action classes
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&) { return Action::OK; }
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    DataFilter dataFilter_; ///< Class used to filter data set(s)
};
#endif
