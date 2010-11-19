#ifndef INC_ACTION_IMAGE_H
#define INC_ACTION_IMAGE_H
// Image
#include "Action.h"

class Image: public Action {
    AtomMask Mask1;
    AtomMask *ComMask;
    bool origin;
    bool center;
    bool ortho;
    bool useMass;
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
