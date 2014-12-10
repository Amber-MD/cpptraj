#ifndef INC_ACTION_CENTER_H
#define INC_ACTION_CENTER_H
/// Action to center coordinates to coord origin or box origin.
#include "Action.h"
class Action_Center: public Action {
  public:
    Action_Center();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Center(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    AtomMask Mask_;
    Frame::CenterMode centerMode_;
    bool useMass_;
    Vec3 refCenter_;
};
#endif
