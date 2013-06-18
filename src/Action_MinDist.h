#ifndef INC_ACTION_MINDIST_H
#define INC_ACTION_MINDIST_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_MinDist
/// Action to calculate the min. distance between atoms in two masks.
class Action_MinDist: public Action {
  public:
    Action_MinDist();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MinDist(); }
    static void Help();
  private:
    ImagedAction image_;
    DataSet* dist_;  ///< Will hold DataSet of minimum distances.
    AtomMask Mask1_;
    AtomMask Mask2_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}
};
#endif  
