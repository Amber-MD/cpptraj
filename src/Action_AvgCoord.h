#ifndef INC_ACTION_AVGCOORD_H
#define INC_ACTION_AVGCOORD_H
/// Class: Action_AvgCoord
/// Action to calculate the overall average of all atomic coords. 
#include "Action.h"
class Action_AvgCoord: public Action {
  public:
    Action_AvgCoord();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_AvgCoord(); }
    static void Help();

    ~Action_AvgCoord();

    void Print() {}
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    bool calcMagnitude_;
    bool useMass_;
    AtomMask Mask_;
    CpptrajFile outfile_;
};
#endif
