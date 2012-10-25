#ifndef INC_ACTION_OUTTRAJ_H
#define INC_ACTION_OUTTRAJ_H
// Action_Outtraj
#include "Action.h"
#include "Trajout.h"
/// Write out a trajectory inside the ActionList
class Action_Outtraj: public Action {
  public:
    Action_Outtraj();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Outtraj(); }
    static void Help();

  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    Trajout outtraj_;
    std::vector<double> Max_;
    std::vector<double> Min_;
    std::vector<DataSet*> Dsets_;
    DataSet* maxmin_;
    Topology* CurrentParm_;
};
#endif
