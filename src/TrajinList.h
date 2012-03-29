#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "CoordFileList.h"
// Class: TrajinList
/// Hold input trajectories
class TrajinList : public CoordFileList {
  public:
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(char*, ArgList *, Topology *);
    /// Set up frames to be processed 
    int SetupFrames();
};
#endif

