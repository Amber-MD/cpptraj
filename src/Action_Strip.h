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

    Topology* newParm_;         ///< Hold stripped topology.
    DataSetList* masterDSL_;    ///< Pointer to master data set list.
    Frame newFrame_;            ///< Hold stripped frame.
    AtomMask M1_;               ///< Atoms to keep.
    ActionTopWriter topWriter_; ///< Used to write stripped topology.
    double charge_;             ///< Total charge to scale remainiing atoms to.
    bool redist_charge_;        ///< If true, rescale charge of remaining atoms to charge_.
};
#endif
