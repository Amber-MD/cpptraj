#ifndef INC_ACTION_DIPOLE_H
#define INC_ACTION_DIPOLE_H
#include "Action.h"
#include "DataSet_Grid.h"
#include "GridAction.h"
class Action_Dipole : public Action, private GridAction {
  public:
    Action_Dipole();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Dipole(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    DataSet_GridFlt* grid_;
    std::vector<double> dipolex_;
    std::vector<double> dipoley_;
    std::vector<double> dipolez_;
    std::string filename_;
    AtomMask mask_;
    double max_;
    Topology* CurrentParm_;
};
#endif
