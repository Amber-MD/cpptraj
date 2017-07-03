#ifndef INC_ACTION_DISTANCE_H
#define INC_ACTION_DISTANCE_H
#include "Action.h"
#include "ImagedAction.h"
/// Action to calculate a distance between atoms in two masks.
class Action_Distance: public Action {
  public:
    Action_Distance();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Distance(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    enum ModeType { NORMAL = 0, REF, PLANE };
    enum PlaneType { YZ = 0, XZ, XY };

    AtomMask Mask1_;     ///< Mask selecting first point
    AtomMask Mask2_;     ///< Mask selecting second point
    ImagedAction image_; ///< Imaging routines.
    DataSet* dist_;      ///< Will hold DataSet of calculated distances.
    ModeType mode_;      ///< Type of distance calculation.
    PlaneType plane_;    ///< If PLANE, which plane to calc distance to.
    bool useMass_;       ///< If true, mass-weight distances.
};
#endif  
