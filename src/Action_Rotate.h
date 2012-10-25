#ifndef INC_ACTION_ROTATE_H
#define INC_ACTION_ROTATE_H
#include "Action.h"
class Action_Rotate : public Action {
  public:
    Action_Rotate();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Rotate(); }
    static void Help();

  private:
    double RotMatrix_[9];
    AtomMask mask_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}
};
#endif
