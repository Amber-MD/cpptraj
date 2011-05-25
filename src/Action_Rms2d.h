#ifndef INC_ACTION_RMS2D_H
#define INC_ACTION_RMS2D_H
// Rms2d
#include "Action.h"

class Rms2d: public Action {
    FrameList ReferenceFrames;
    bool nofit;
    bool useMass;
    AtomMask RefMask;
    AtomMask FrameMask;
    char *rmsdFile;
    DataSetList RmsData;

    void progressBar(int, int);
  public:
    Rms2d();
    ~Rms2d();

    int init();
    int setup();
    int action();
    void print();
};
#endif
