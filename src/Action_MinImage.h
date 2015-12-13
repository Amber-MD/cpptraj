#ifndef INC_ACTION_MINIMAGE_H
#define INC_ACTION_MINIMAGE_H
#include "Action.h"
#include "ImagedAction.h"
//#incl ude "PDBfile.h" // DEBUG
/// Action to calculate minimum non-self distance between atoms in two masks.
class Action_MinImage: public Action {
  public:
    Action_MinImage();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_MinImage(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    double MinNonSelfDist2(Vec3 const&, Vec3 const&);

    ImagedAction image_;
    Matrix_3x3 ucell_, recip_;
    DataSet* dist_;      ///< Will hold DataSet of calculated distances.
    DataSet* atom1_;
    DataSet* atom2_;
    bool useMass_;       ///< If true, mass-weight distances.
    bool calcUsingMask_; ///< If true use center of masks
    AtomMask Mask1_;
    AtomMask Mask2_;
    std::vector<double> minDist_;
    std::vector<int> minAtom1_;
    std::vector<int> minAtom2_;
    //PDBfile pdbout_; // DEBUG
};
#endif
