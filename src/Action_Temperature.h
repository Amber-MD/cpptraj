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

    /// How temperature will be determined.
    enum ModeType {
      FROM_FRAME = 0,  ///< From incoming trajectory frame.
      CALC_ONLY,       ///< From incoming velocities.
      CALC_AND_MODIFY  ///< From velocities, update temperature in frame.
    };

    DataSet* Tdata_;       ///< Hold temperature data
    AtomMask Mask_;
    Constraints cons_;
    CoordinateInfo cInfo_; ///< In case we are modifying Frame temperature
    ModeType mode_;        ///< How we are calculating temperature
    int dof_offset_;       ///< # of global d.o.f. removed
    bool removeTrans_;     ///< True if accounting for removed translational d.o.f.
    bool removeRot_;       ///< True if accounting for removed rotational d.o.f.
};
#endif
