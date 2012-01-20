#ifndef INC_ACTION_DISTRMSD_H
#define INC_ACTION_DISTRMSD_H
#include "Action.h"
#include "TrajectoryFile.h"
// Class: DistRmsd
/// Action to calculate the distance RMSD between frame and a reference frame.
class DistRmsd: public Action {
    DataSet *drmsd;
    AtomMask RefMask, TgtMask;    ///< Target and reference masks.
    Frame RefFrame;               ///< Hold reference frame coords
    Frame SelectedRef;            ///< Hold only ref coods selected by RefMask 
    Frame SelectedTgt;            ///< Hold only target coords selected by TgtMask
    TrajectoryFile *RefTraj;      ///< Reference trajectory, each frame used in turn
    AmberParm *RefParm;           ///< Reference frame Parm
    bool first;                   ///< If true, use first frame read in as reference.

    int SetRefMask();

  public:
    DistRmsd();
    ~DistRmsd();

    int init();
    int setup();
    int action();
};
#endif
