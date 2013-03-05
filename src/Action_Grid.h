#ifndef INC_ACTION_GRID_H
#define INC_ACTION_GRID_H
#include "Action.h"
#include "Grid.h"
class Action_Grid : public Action {
  public:
    Action_Grid();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Grid(); }
    static void Help();


    void Print();
  private:
    double max_;
    double madura_;
    double smooth_;
    bool invert_;
    bool dxform_;
    AtomMask mask_;
    std::string filename_;
    std::string pdbname_;
    Grid grid_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
};
#endif
