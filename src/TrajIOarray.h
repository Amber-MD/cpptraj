#ifndef INC_TRAJIOARRAY_H
#define INC_TRAJIOARRAY_H
#include "TrajectoryIO.h"
#include "TrajFrameCounter.h"
/// Hold an array of TrajectoryIO classes.
class TrajIOarray {
  public:
    TrajIOarray() : debug_(0) {}
    ~TrajIOarray();
    /// Clear all IO classes and names.
    void ClearIOarray();
    /// Set up replica file names from arguments or by searching.
    int SetupReplicaFilenames(std::string const&, ArgList&);
    /// Given lowest replica name, search for all other replicas.
    int SearchForReplicas(std::string const&);
    /// Add the given file name and names from comma-separated list. 
    int AddReplicasFromArgs(std::string const&, std::string const&);
    /// Setup TrajectoryIO classes for all file names.
    int SetupIOarray(ArgList&, TrajFrameCounter&, CoordinateInfo&, Topology*);
    /// Print file names and IO info to STDOUT
    void PrintIOinfo() const;
    /// Iterator over trajectories
    typedef std::vector<TrajectoryIO*>::const_iterator const_iterator;
    /// Iterator to beginning of array.
    const_iterator begin() const { return IOarray_.begin(); }
    /// Iterator to end of array
    const_iterator end()   const { return IOarray_.end();   }
    /// \return specified TrajectoryIO pointer
    TrajectoryIO* operator[](unsigned int u) const { return IOarray_[u]; }
    /// \return number of set-up trajectories
    size_t size()          const { return IOarray_.size();  }
    /// \return specified file name.
    std::string const& F_name(unsigned int u) const { return replica_filenames_[u]; }
    const char* f_name(unsigned int u) const { return replica_filenames_[u].c_str(); }
    /// 'remdout' deprecated error message.
    static const char* DEPRECATED_remdout;
  private:
    typedef std::vector<TrajectoryIO*> IOarrayType;
    typedef std::vector<std::string> NameListType;

    IOarrayType IOarray_;           ///< Input replica trajectories.
    NameListType replica_filenames_; ///< Replica traj file names.
    int debug_;
};
#endif
