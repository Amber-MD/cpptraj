#ifndef INC_ACTION_CHANNEL_H
#define INC_ACTION_CHANNEL_H
#include "Action.h"
/// Experimental action for calculating solvent channels.
class Action_Channel : public Action {
  public:
    Action_Channel();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Channel(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet* grid_;
    AtomMask soluteMask_;
    AtomMask solventMask_;
    Vec3 dxyz_;
    std::vector<double> radii_;
};
#endif
