#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
#include <string>
#include <vector>
// Forward declarations
class Trajout_Single;
class Topology;
class ArgList;
class DataSetList;
class CoordinateInfo;
class Frame;
/// Hold output trajectories for a run.
class TrajoutList {
  public:
    TrajoutList();
    ~TrajoutList();
    void SetDebug(int);
    void Clear();
    /// FIXME legacy function to maintain pytraj compatibility
    int AddTrajout(std::string const&, ArgList const&, Topology*);
    /// Add output trajectory to the list and associate with given topology.
    int AddTrajout(std::string const&, ArgList const&, DataSetList const&, Topology*);
    /// Set up trajectories for given topology.
    int SetupTrajout(Topology*, CoordinateInfo const&, int);
    /// Write frame to normal output trajectories.
    int WriteTrajout(int, Frame const&);
    /// Call end for all trajectories
    void CloseTrajout();
    /// List output trajectories.
    void List(std::vector<int> const&) const;
    /// \return true if no args/trajectories present.
    bool Empty() const { return trajout_.empty(); }
#   ifdef MPI
    int ParallelSetupTrajout(Topology*, CoordinateInfo const&, int, Parallel::Comm const&);
#   endif
  private:
    void ListActive() const;

    int debug_;
    typedef std::vector<Topology*> TopArray;
    typedef std::vector<Trajout_Single*> ListType;
    ListType trajout_;       ///< Hold output trajectories.
    ListType active_;        ///< Hold only active output trajectories.
    TopArray trajoutTops_;   ///< Array of associated topology files
    std::vector<bool> open_; ///< True if corresponding trajectory is open
};
#endif
