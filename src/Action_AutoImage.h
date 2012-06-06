#ifndef INC_ACTION_AUTOIMAGE_H
#define INC_ACTION_AUTOIMAGE_H
#include "Action.h"
class Action_AutoImage : public Action {
  public:
    Action_AutoImage();

  private:
    int init();
    int setup();
    int action();

    // Reference vars
    enum RefModeType {UNKNOWN_REF=0, FIRST, REF};
    RefModeType refmode_;
    Frame RefFrame_;
    Frame SeledctedRef_;
    
    AtomMask anchorMask_; ///< Used to center anchor region
    std::string anchor_;  ///< Mask expression for anchor region
    std::string fixed_;   ///< Mask expression for fixed region
    std::string mobile_;  ///< Mask expression for mobile region

    bool origin_;
    bool ortho_;
    bool center_;
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic_;

    std::vector<int> anchorList_;
    std::vector<int> fixedList_;
    std::vector<int> mobileList_;

    std::vector<int> SetupAtomRanges(std::string const&);
};
#endif
