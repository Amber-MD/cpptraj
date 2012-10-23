#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "TrajectoryFile.h"
// Class: TrajinList
/// Hold input trajectories
class TrajinList {
  public:
    TrajinList();
    ~TrajinList();
    static void Help();
    void SetDebug(int dIn) { debug_ = dIn; }
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(ArgList*, Topology*);
    void Begin();
    TrajectoryFile *NextTraj();
    void List();
    int MaxFrames() { return maxframes_; }
  private:
    std::vector<TrajectoryFile*> trajin_;
    int debug_;
    int maxframes_;
    std::vector<TrajectoryFile*>::iterator currentTraj_;
};
#endif

