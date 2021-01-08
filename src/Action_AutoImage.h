#ifndef INC_ACTION_AUTOIMAGE_H
#define INC_ACTION_AUTOIMAGE_H
#include "Action.h"
namespace Image {
  class List_Unit;
}
/// Perform imaging, attempting to keep solute in 1 configuration
class Action_AutoImage : public Action {
  public:
    Action_AutoImage();
    ~Action_AutoImage();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AutoImage(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    AtomMask anchorMask_; ///< Used to center anchor region.
    std::string anchor_;  ///< Mask expression for anchor region.
    std::string fixed_;   ///< Mask expression for fixed region.
    std::string mobile_;  ///< Mask expression for mobile region.
    int debug_;
    bool origin_;         ///< If true imaging occurs w.r.t. coordinate origin.
    bool usecom_;         ///< If true imaging of mobile region uses molecule center.
    bool truncoct_;       ///< If true image into truncated octahedron shape.
    bool useMass_;        ///< If true use center of mass
    bool movingAnchor_;   ///< If true anchor position set to previous fixed molecule
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic_; ///< Determine whether triclinic code should be used.

    Image::List_Unit* fixedList_;  ///< Contain atom indices for fixed elements
    Image::List_Unit* mobileList_; ///< Contain atom indices for mobile elements.
};
#endif
