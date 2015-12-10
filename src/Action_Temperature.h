#ifndef INC_ACTION_TEMPERATURE_H
#define INC_ACTION_TEMPERATURE_H
#include "Action.h"
/// Calculate the temperature of parts of a system.
class Action_Temperature : public Action {
  public:
    Action_Temperature();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Temperature(); }
    void Help() const;
  private:
    enum ShakeType {OFF = 0, BONDS_TO_H, ALL_BONDS};
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet* Tdata_;
    bool getTempFromFrame_;
    AtomMask Mask_;
    ShakeType shakeType_;
    int degrees_of_freedom_;
};
#endif
