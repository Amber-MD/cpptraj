#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "TrajectoryFile.h"
// Class: TrajinList
/// Hold input trajectories
class TrajinList {
  public:
    TrajinList();
    ~TrajinList();
    void SetDebug(int dIn) { debug_ = dIn; }
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(ArgList*, Topology*);
    /// Set up frames to be processed 
    int SetupFrames();

    void Begin();
    TrajectoryFile *NextTraj();
  private:
    std::vector<TrajectoryFile*> trajin_;
    int debug_;
    std::vector<TrajectoryFile*>::iterator currentTraj_;
};
#endif

