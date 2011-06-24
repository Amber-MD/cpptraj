#ifndef INC_ACTION_RMS2D_H
#define INC_ACTION_RMS2D_H
// Rms2d
#include "Action.h"
#include "TrajectoryFile.h"

class Rms2d: public Action {
    FrameList ReferenceFrames; // Hold frames from all trajin stmts
    bool nofit;                // Do not perform rms fitting
    bool useMass;              // Perform mass-weighted rmsd
    AtomMask RefMask;          // Reference atom mask
    AtomMask FrameMask;        // Target atom mask
    char *rmsdFile;            // Output filename
    DataSetList RmsData;       // 1 data set for each ref frame to each tgt frame
    TrajectoryFile *RefTraj;   // Reference trajectory, each frame used in turn
    AmberParm *RefParm;        // Reference trajectory Parm

  public:
    Rms2d();
    ~Rms2d();

    int init();
    int setup();
    int action();
    void print();
};
#endif
