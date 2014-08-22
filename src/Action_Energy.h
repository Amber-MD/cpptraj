#ifndef INC_ACTION_ENERGY_H
#define INC_ACTION_ENERGY_H
#include "Action.h"
#include "Energy.h"
/// Calculate energy 
class Action_Energy: public Action {
  public:
    Action_Energy();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Energy(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();

    enum Etype { BOND, ANGLE, DIHEDRAL, V14, Q14, VDW, ELEC, TOTAL};
    std::vector<DataSet*> Energy_;
    Topology* currentParm_;
    AtomMask Mask1_;
    AtomMask Imask_;
    Energy_Amber ENE_;
};
#endif
