#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
#include "Trajout_Single.h"
#include "EnsembleOut.h"
/// Hold output ensembles.
class EnsembleOutList {
  public:
    EnsembleOutList() {}
    ~EnsembleOutList() { Clear(); } 
    void Clear();
    int AddEnsembleOut(std::string const&, ArgList const&, Topology*,
                       int, TrajectoryFile::TrajFormatType);
    int SetupEnsembleOut(Topology*, CoordinateInfo const&, int); //TODO const, Topology array?
    int WriteEnsembleOut(int, FramePtrArray const&);
    void CloseEnsembleOut();
    void List(std::vector<int> const&) const;
  private:
    typedef std::vector<Topology*> TopArray;
    typedef std::vector<EnsembleOut*> EnsArray;
    EnsArray ensout_;
    EnsArray active_;
    TopArray ensTops_;
    std::vector<bool> open_;
};
// =============================================================================
/// Hold output trajectories for a run.
/** When a 'trajout' or equivalent command is given, the trajectory will
  * be added to TrajoutList as a Trajout_Single class by default. This is
  * done so that some error checking of arguments etc can occur, even if
  * eventually the list will be used in ensemble output.
  */
class TrajoutList {
  public:
    TrajoutList() : debug_(0) {}
    ~TrajoutList() { Clear(); }
    void SetDebug(int);
    void Clear();
    /// Add output trajectory to the list and associate with given topology.
    int AddTrajout(std::string const&, ArgList const&, Topology*);
    /// \return Array with current output trajectories converted to ensemble output trajectories.
    int MakeEnsembleTrajout(EnsembleOutList&, int) const;
    /// Set up trajectories for given topology.
    int SetupTrajout(Topology*, CoordinateInfo const&, int);
    /// Write frame to normal output trajectories.
    int WriteTrajout(int, Frame const&);
    /// Call end for all trajectories
    void CloseTrajout();
    /// List output trajectories.
    void List(std::vector<int> const&) const;
    /// \return true if no args/trajectories present.
    bool Empty()     const { return trajout_.empty();     }
#   ifdef MPI
    int ParallelSetupTrajout(Topology*, CoordinateInfo const&, int, Parallel::Comm const&);
#   endif
  private:
    int debug_;
    typedef std::vector<Trajout_Single*> ListType;
    typedef std::vector<ArgList> ArgsArray;
    typedef std::vector<Topology*> TopArray;
    typedef std::vector<std::string> Sarray;
    ListType trajout_; ///< Hold output trajectories.
    ListType active_;  ///< Hold only active output trajectories.
    ArgsArray trajoutArgs_; ///< Array of trajout args for potentially setting up ensemble.
    TopArray  trajoutTops_; ///< Array of associated topologies.
    Sarray trajoutNames_;   ///< Array of trajout file names.
    std::vector<bool> open_;
};
#endif
