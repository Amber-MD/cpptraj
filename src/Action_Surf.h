#ifndef INC_ACTION_SURF_H
#define INC_ACTION_SURF_H
#include "Action.h"
// Class: Action_Surf
/// Calculate LCPO surface area.
/** LCPO method from:
  * -  J. Weiser, P.S. Shenkin, and W.C. Still,
  *    "Approximate atomic surfaces from linear combinations of pairwise
  *    overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
  */
class Action_Surf: public Action {
  public:
    Action_Surf();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Surf(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}
    /// Assign LCPO vdW radius and parameters 1-4.
//    void AssignLCPO(SurfInfo*, double, double, double, double, double);
//    void SetAtomLCPO(Topology const&,int, SurfInfo*);

    /// Contain LCPO parameters
    class SurfInfo;

    DataSet* surf_;        ///< Hold LCPO surface area data
    AtomMask Mask1_;       ///< Atoms to calculate SA for
    AtomMask SoluteAtoms_; ///< All solute atoms

};
/// Hold LCPO parameters
class Action_Surf::SurfInfo {
  public:
    SurfInfo() : vdw_(0.0), P1_(0.0), P2_(0.0), P3_(0.0), P4_(0.0) {}
//    SurfInfo(double v, double p1, double p2, double p3, double p4) :
//      vdw_(v), P1_(p1), P2_(p2), P3_(p3), P4_(p4) {}
    /// CONSTRUCTOR - Set LCPO parameters from atom element/type/#bonds
    SurfInfo(Atom const&);
    double VDW() const { return vdw_; }
    double P1()  const { return P1_; }
    double P2()  const { return P2_; }
    double P3()  const { return P3_; }
    double P4()  const { return P4_; }
  private:
    double vdw_;
    double P1_;
    double P2_;
    double P3_;
    double P4_;
};
#endif
