#ifndef INC_ACTION_FIXATOMORDER_H
#define INC_ACTION_FIXATOMORDER_H
#include "Action.h"
#include "ActionTopWriter.h"
#include "CharMask.h"
// Class: Action_FixAtomOrder
/// Fix atom ordering in parm where atoms in mols are not sequential. 
class Action_FixAtomOrder: public Action {
  public:
    Action_FixAtomOrder();
    ~Action_FixAtomOrder();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_FixAtomOrder(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    void VisitAtom(int,int,Topology const&);
    Action::RetType PdbOrder(ActionSetup&);
    Action::RetType FixMolecules(ActionSetup&);

    enum ModeType {FIX_MOLECULES = 0, PDB_ORDER};

    int debug_;
    typedef std::vector<int> MapType;
    MapType atomMap_;           ///< Map original atoms to new atoms.
    Topology* newParm_;         ///< Re-ordered topology
    Frame newFrame_;            ///< Re-ordered frame
    ActionTopWriter topWriter_; ///< Use to write re-ordered Topology
    ModeType mode_;             ///< How atoms should be reordered
    CharMask hetatm_;           ///< If PDB reorder, atoms that should be HETATM
};
#endif
