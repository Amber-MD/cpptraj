#ifndef INC_ACTION_DIHEDRAL_H
#define INC_ACTION_DIHEDRAL_H
#include "Action.h"
/// Calculate dihedral in a Frame
class Action_Dihedral: public Action {
  public:
    Action_Dihedral();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Dihedral(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
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
