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

    DataSet* dist_;  ///< Will hold DataSet of calculated distances.
    bool useMass_;   ///< If true, mass-weight distances.
    AtomMask Mask1_;
    AtomMask Mask2_;
    ImagedAction image_; ///< Imaging routines.
};
#endif  
