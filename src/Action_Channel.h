#ifndef INC_ACTION_CHANNEL_H
#define INC_ACTION_CHANNEL_H
#include "Action.h"
#include "Grid.h"
class Action_Channel : public Action {
  public:
    Action_Channel();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Channel(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    Grid solvent_;
    Grid solute_;
    AtomMask solventMask_;
    AtomMask soluteMask_;
};
#endif
