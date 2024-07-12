#ifndef INC_ACTION_ADDATOM_H
#define INC_ACTION_ADDATOM_H
#include "Action.h"
#include "ActionTopWriter.h"
/// Add an atom to Topology/coordinates 
class Action_AddAtom : public Action {
  public:
    /// CONSTRUCTOR
    Action_AddAtom();
    /// DESTRUCTOR
    ~Action_AddAtom();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AddAtom(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    Topology* newParm_;         ///< Hold topology with added atom
    Frame newFrame_;            ///< Hold coords with added atom
    ActionTopWriter topWriter_; ///< Used to write new topology.
    Atom newAtom_;              ///< New atom
    NameType residueName_;      ///< Residue name new atom will belong to.
    Vec3 xyz_;                  ///< New atom position
};
#endif
