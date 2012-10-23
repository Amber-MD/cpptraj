#ifndef INC_ACTION_CENTER_H
#define INC_ACTION_CENTER_H
/// Class: Action_Center
/// Action to center coordinates to coord origin or box origin.
#include "Action.h"
class Action_Center: public Action {
  public:
    Action_Center();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Center(); }
    static void Help();

    void Print() {}
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    AtomMask Mask_;
    bool origin_;
    bool useMass_;
};
#endif
