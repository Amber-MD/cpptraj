#ifndef INC_ACTION_AREAPERMOL_H
#define INC_ACTION_AREAPERMOL_H
#include "Action.h"
/// Calculate area per molecule for given molecule/dimensions. 
class Action_AreaPerMol: public Action {
  public:
    Action_AreaPerMol();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AreaPerMol(); }
    void Help() const;
  private:
    enum AreaType { XY, XZ, YZ };

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    DataSet *area_per_mol_;
    double Nmols_;
    double Nlayers_;
    AreaType areaType_;
    CharMask Mask1_;
};
#endif
