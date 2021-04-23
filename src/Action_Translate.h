#ifndef INC_ACTION_TRANSLATE_H
#define INC_ACTION_TRANSLATE_H
#include "Action.h"
class Action_Translate : public Action {
  public:
    Action_Translate();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Translate(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Vec3 Trans_;    ///< Hold X/Y/Z translation, or point to translate to.
    AtomMask mask_; ///< Mask of atoms to translate.
    bool toPoint_;  ///< If true, translate to a point instead of by deltas.
    bool useMass_;  ///< If true, to-point translation weighted by mass.
};
#endif
