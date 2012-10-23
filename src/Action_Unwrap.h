#ifndef INC_ACTION_UNWRAP_H
#define INC_ACTION_UNWRAP_H
#include "Action.h"
class Action_Unwrap : public Action {
  public:
    Action_Unwrap();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Unwrap(); }
    static void Help();

    void Print() {}
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    AtomMask mask_;
    Frame RefFrame_;
    Topology* RefParm_;
    bool orthogonal_;
};
#endif
