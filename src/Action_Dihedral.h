#ifndef INC_ACTION_DIHEDRAL_H
#define INC_ACTION_DIHEDRAL_H
#include "Action.h"
/// Calculate dihedral in a Frame
class Action_Dihedral: public Action {
  public:
    Action_Dihedral();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Dihedral(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet* dih_;
    double minTorsion_; ///< Values less than this will be shifted +360.0
    bool useMass_;
    AtomMask M1_;
    AtomMask M2_;
    AtomMask M3_;
    AtomMask M4_;
};
#endif  
