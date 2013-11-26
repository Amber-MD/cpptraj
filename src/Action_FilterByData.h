#ifndef INC_ACTION_FILTERBYDATA_H
#define INC_ACTION_FILTERBYDATA_H
#include "Action.h"
#include "Array1D.h"
/// Filter out frames by DataSet
class Action_FilterByData : public Action {
  public:
    Action_FilterByData() : maxmin_(0) {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_FilterByData(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**) { return Action::OK; }
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    std::vector<double> Max_;
    std::vector<double> Min_;
    Array1D Dsets_;
    DataSet* maxmin_;
};
#endif
