#ifndef INC_ACTION_TEMPERATURE_H
#define INC_ACTION_TEMPERATURE_H
#include "Action.h"
#include "Constraints.h"
/// Calculate the temperature of parts of a system.
class Action_Temperature : public Action {
  public:
    Action_Temperature();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Temperature(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet* Tdata_;
    bool getTempFromFrame_;
    AtomMask Mask_;
    Constraints cons_;
    bool removeTrans_; ///< True if accounting for removed translational d.o.f.
    bool removeRot_;   ///< True if accounting for removed rotational d.o.f.
    int dof_offset_;   ///< # of global d.o.f. removed
};
#endif
