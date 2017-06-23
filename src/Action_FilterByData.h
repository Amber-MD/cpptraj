#ifndef INC_ACTION_FILTERBYDATA_H
#define INC_ACTION_FILTERBYDATA_H
#include "Action.h"
#include "Array1D.h"
/// Filter out frames by DataSet
class Action_FilterByData : public Action {
  public:
    Action_FilterByData() : maxmin_(0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_FilterByData(); }
    void Help() const;
    /// For running as a separate command.
    size_t DetermineFrames() const;
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType DoAction(int, ActionFrame&);
  private:
    Action::RetType Setup(ActionSetup&) { return Action::OK; }
    void Print();

    typedef std::vector<DataSet*> DSarray;

    std::vector<double> Max_;
    std::vector<double> Min_;
    Array1D Dsets_;
    DSarray outsets_;
    DataSet* maxmin_;
    int Npassed_;
    int Nfiltered_;
    bool multi_;
};
#endif
