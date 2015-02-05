#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
#include "Trajout.h"
#include "TopologyList.h"
// Class: TrajoutList
/// Hold trajectories for output
class TrajoutList {
  public:
    TrajoutList();
    ~TrajoutList();
    void Clear();
    void SetDebug(int);
    /// Place current trajout in given list as ensemble trajectory output.
    int MakeEnsembleTrajout(TopologyList const&, TrajoutList&);
    /// Add trajectory to list for single trajectory output
    int AddTrajout(ArgList const&, TopologyList const&);
    /// Write frame array to ensemble output trajectories.
    int WriteEnsembleOut(int, FramePtrArray const&);
    /// Set up trajectories for given topology.
    int SetupTrajout(Topology*);
    /// Write frame to normal output trajectories.
    int WriteTrajout(int, Frame const&);
    /// Call end for all trajectories
    void CloseTrajout();
    /// List output trajectories.
    void List() const;
    /// \return true if no args/trajectories present.
    bool Empty()     const { return trajout_.empty();     }
  private:
    int debug_;
    typedef std::vector<Trajout*> ListType;
    ListType trajout_; ///< Hold actual output trajectories.
    ListType active_;  ///< Hold only active output trajectories.
    typedef std::vector<ArgList> ArgsArray;
    ArgsArray trajoutArgs_; ///< Array of trajout args for setting up trajouts.
};
#endif
