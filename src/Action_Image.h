#ifndef INC_ACTION_IMAGE_H
#define INC_ACTION_IMAGE_H
/// Action to wrap coordinates back into primary box
#include "Action.h"
#include "ImageTypes.h"
class Action_Image: public Action {
  public:
    Action_Image();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Image(); }
    void Help() const;
    ~Action_Image();
  private:
    Image::Mode imageMode_;
    /// Mask expression selecting atoms to image.
    std::string maskExpression_;
    /// If defined, image w.r.t. the center of atoms in ComMask.
    AtomMask *ComMask_;
    /// Offsets
    Vec3 offset_;
    /// If true image w.r.t. coordinate origin, otherwise box center
    bool origin_;
    /// If true molecules will be imaged w.r.t. their center, otherwise first atom will be used
    bool center_;
    /// True if orthorhombic cell, false otherwise.
    bool ortho_;
    bool useMass_;
    bool truncoct_;
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic_;
    int debug_;
    /// Vector containing atom ranges to be imaged (first to last)
    Image::PairType imageList_; 

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}
};
#endif
