#ifndef INC_TRAJIOARRAY_H
#define INC_TRAJIOARRAY_H
#include "TrajectoryIO.h"
#include "TrajFrameCounter.h"
#include "FileName.h"
/// Hold an array of TrajectoryIO classes.
class TrajIOarray {
  public:
    TrajIOarray() : debug_(0) {}
    ~TrajIOarray();
    /// Clear all IO classes and names.
    void ClearIOarray();
    /// Set up replica file names from arguments or by searching.
    int SetupReplicaFilenames(FileName const&, ArgList&);
    /// Setup TrajectoryIO classes for all file names.
    int SetupIOarray(ArgList&, TrajFrameCounter&, CoordinateInfo&, Topology*);
#   ifdef MPI
    /// Set up replica filenames such that each ensemble rank checks 1 file
    int SetupReplicaFilenames(FileName const&, ArgList&, Parallel::Comm const&,
                              Parallel::Comm const&);
    /// Each ensemble rank sets up TrajectoryIO class only for member it will process.
    int SetupIOarray(ArgList& argIn, TrajFrameCounter& counter,
                     CoordinateInfo& cInfo, Topology* trajParm,
                     Parallel::Comm const&, Parallel::Comm const&);
#   endif
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
    /// Set debug level
    void SetDebug(int d) { debug_ = d; }
  private:
    /// Add the given file name and names from comma-separated list. 
    int AddReplicasFromArgs(FileName const&, std::string const&);
#   ifdef MPI
    /// Each rank searches for replica from given name + comma-separated list
    int AddReplicasFromArgs(FileName const&, std::string const&, Parallel::Comm const&,
                            Parallel::Comm const&);
#   endif
    typedef std::vector<TrajectoryIO*> IOarrayType;
    IOarrayType IOarray_;               ///< Input replica trajectories.
    File::NameArray replica_filenames_; ///< Replica traj file names.
    int debug_;
};
#endif
