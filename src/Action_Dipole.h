#ifndef INC_ACTION_DIPOLE_H
#define INC_ACTION_DIPOLE_H
#include "Action.h"
#include "Grid.h"
class Action_Dipole : public Action {
  public:
    Action_Dipole();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Dipole(); }
    static void Help();

    void Print();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    Grid grid_;
    std::vector<double> dipolex_;
    std::vector<double> dipoley_;
    std::vector<double> dipolez_;
    std::string filename_;
    AtomMask mask_;
    double max_;
    Topology* CurrentParm_;
};
#endif
