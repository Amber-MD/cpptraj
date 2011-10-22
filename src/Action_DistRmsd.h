#ifndef INC_ACTION_DISTRMSD_H
#define INC_ACTION_DISTRMSD_H
/// Class: DistRmsd
/// Action to calculate the distance RMSD between frame and a reference frame.
#include "Action.h"
#include "TrajectoryFile.h"
class DistRmsd: public Action {
    DataSet *drmsd;
    AtomMask RefMask, TgtMask;          // Target and reference masks.
    bool first;                         // Action options
    Frame RefFrame;                    // Hold reference frame coords
    Frame SelectedRef;                 // Hold only ref coods selected by RefMask 
    Frame SelectedTgt;                 // Hold only target coords selected by TgtMask
    TrajectoryFile *RefTraj;            // Reference trajectory, each frame used in turn
    AmberParm *RefParm;                 // Reference frame Parm

    int SetRefMask();

  public:
    DistRmsd();
    ~DistRmsd();

    int init();
    int setup();
    int action();
};
#endif
