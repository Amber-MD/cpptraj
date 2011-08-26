#ifndef INC_ACTION_IMAGE_H
#define INC_ACTION_IMAGE_H
/// Class: Image
/// Action to wrap coordinates back into primary box
#include "Action.h"
class Image: public Action {
    AtomMask Mask1;
    AtomMask *ComMask;
    bool origin;
    bool center;
    bool ortho;
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic;

  public:
    Image();
    ~Image();

    int init();
    int setup();
    int action();
};
#endif
