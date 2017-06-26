#ifndef INC_ACTION_UNIMAGE_H
#define INC_ACTION_UNIMAGE_H
#include "Action.h"
#include "ImagedAction.h"
#include "ImageTypes.h"
/// <Enter description of Action_Unimage here>
class Action_Unimage : public Action {
  public:
    Action_Unimage() : natoms_(0), useCenter_(false) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Unimage(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    typedef Image::PairType Iarray;

    ImagedAction image_ ; ///< Imaging routines
    Frame previous_; ///< Hold previous frames selected coords
    Iarray atomPairs_; ///< Hold start and end atom indices for each element in current to unwrap
    Vec3 boxCenter_; ///< Box center for current frame
    Vec3 boxTrans_; ///< Translation vector for current element
    Matrix_3x3 ucell_; ///< Unit cell matrix for current frame
    Matrix_3x3 recip_; ///< Recip (frac) matrix for current frame
    std::string maskExp_; ///< Mask expression
    int natoms_; ///< Total # atoms being unwrapped
    Image::Mode imageMode_; ///< Unwrap mode
    bool useCenter_; ///< If true use c.o.m. of elements, otherwise first atom

};
#endif
