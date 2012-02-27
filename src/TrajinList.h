#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "CoordFileList.h"
// Class: TrajinList
/// Hold input trajectories
class TrajinList : public CoordFileList {
  public:
    TrajinList();
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(char*, ArgList *, AmberParm *);
    /// Set up frames to be processed 
    int SetupFrames();
};
#endif

