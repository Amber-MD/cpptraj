#ifndef INC_ACTION_SETVELOCITY_H
#define INC_ACTION_SETVELOCITY_H
#include "Action.h"
#include "Random.h"
#include "Constraints.h"
/// Set velocities for selected atoms in a system. 
class Action_SetVelocity : public Action {
  public:
    Action_SetVelocity();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_SetVelocity(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    typedef std::vector<double> Darray;

    AtomMask Mask_;    ///< Atoms to set velocities of
    Darray SD_;        ///< Hold sqrt(kB*(1/mass)) for each atom
    double tempi_;     ///< Temperature to generate velocity distribution at.
    Constraints cons_; ///< Hold constraint info 
    Random_Number RN_;
    CoordinateInfo cInfo_;
    Frame newFrame_;
    bool zeroMomentum_;
};
#endif
