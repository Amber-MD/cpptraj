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
    int SetupReplicaFilenames(FileName const&, ArgList&);
#   ifdef MPI
    int SetupReplicaFilenames(FileName const&, ArgList&, Parallel::Comm const&,
                              Parallel::Comm const&);
#   endif
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
    FileName const& F_name(unsigned int u) const { return replica_filenames_[u];        }
    const char* f_name(unsigned int u)     const { return replica_filenames_[u].full(); }
    /// 'remdout' deprecated error message.
    static const char* DEPRECATED_remdout;
  private:
    /// Used for replica filename searching.
    class RepName;
    /// Split given filename into components for searching for other replicas.
    int GetRepPrefixAndExt(FileName const&, std::string&, std::string&, std::string&,
                           int&, int&) const;
    /// Given lowest replica name, search for all other replicas.
    int SearchForReplicas(FileName const&);
    /// Add the given file name and names from comma-separated list. 
    int AddReplicasFromArgs(FileName const&, std::string const&);
#   ifdef MPI
    int SearchForReplicas(FileName const&, Parallel::Comm const&);
    int AddReplicasFromArgs(FileName const&, std::string const&, Parallel::Comm const&);
#   endif
    typedef std::vector<TrajectoryIO*> IOarrayType;
    IOarrayType IOarray_;               ///< Input replica trajectories.
    File::NameArray replica_filenames_; ///< Replica traj file names.
    int debug_;
};
// ----- PRIVATE CLASSES -------------------------------------------------------
/** Given lowest replica traj filename, split into components for search. */
class TrajIOarray::RepName {
  public:
    RepName() : ExtWidth_(0), lowestRepnum_(-1) {}
    RepName(FileName const&, int);
    bool Error() const { return Prefix_.empty(); }
    /// \return Replica file name for given replica number.
    FileName RepFilename(int) const;
    /// \return integer value of lowest replica traj filename extension.
    int LowestRepNum() const { return lowestRepnum_; }
  private:
    std::string Prefix_;      ///< File name up to the numerical extension.
    std::string ReplicaExt_;  ///< Numerical extension.
    std::string CompressExt_; ///< Optional compression extension after numerical extension.
    int ExtWidth_;            ///< Width of the numerical extension. TODO remove
    int lowestRepnum_;        ///< Integer value of numerical extension.
};
#endif
