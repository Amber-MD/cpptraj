#ifndef INC_ACTION_STRIP_H
#define INC_ACTION_STRIP_H
#include "Action.h"
#include "ActionTopWriter.h"
/// Used to remove atoms from the state.
class Action_Strip: public Action {
  public:
    Action_Strip();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Strip(); }
    void Help() const;
    ~Action_Strip();
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Topology* newParm_;
    CoordinateInfo* newCinfo_;
    DataSetList* masterDSL_;
    Frame newFrame_;
    AtomMask M1_;
    bool removeBoxInfo_;
    ActionTopWriter topWriter_;
};
#endif
