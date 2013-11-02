#ifndef INC_ACTION_MAXDIST_H
#define INC_ACTION_MAXDIST_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_MaxDist
/// Action to calculate the max distance between atoms in 1 or 2 masks.
class Action_MaxDist: public Action {
  public:
    Action_MaxDist();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MaxDist(); }
    static void Help();
  private:
    ImagedAction image_;
    DataSet* dist_;  ///< Will hold DataSet of max distances.
    AtomMask Mask1_;
    AtomMask Mask2_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}
};
#endif  
