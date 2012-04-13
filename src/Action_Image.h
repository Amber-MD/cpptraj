#ifndef INC_ACTION_IMAGE_H
#define INC_ACTION_IMAGE_H
// Class: Image
/// Action to wrap coordinates back into primary box
#include "Action.h"
class Image: public Action {
  public:
    Image();
    ~Image();
  private:
    /// Only atoms in Mask1 will be imaged
    AtomMask Mask1;
    /// If defined, image w.r.t. the center of atoms in ComMask.
    AtomMask *ComMask;
    /// If true image w.r.t. coordinate origin, otherwise box center
    bool origin;
    /// If true molecules will be imaged w.r.t. their center, otherwise first atom will be used
    bool center;
    /// True if orthorhombic cell, false otherwise.
    bool ortho;
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic;
    /// Vector containing atom ranges to be imaged (first to last)
    std::vector<int> imageList; 


    int init();
    int setup();
    int action();
};
#endif
