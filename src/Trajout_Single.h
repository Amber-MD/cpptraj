#ifndef INC_TRAJOUT_SINGLE_H
#define INC_TRAJOUT_SINGLE_H
#include "OutputTrajCommon.h"
// Forward declarations
class Frame;
class DataSetList;
/// Write out 1 frame at a time to a single file.
/** Note that unlike Trajin, there is really no point in having
  * a single frame written to multiple files since this is handled
  * by TrajoutList. Therefore this class doesnt need to inherit
  * from a base class currently. The Trajout_Single name is used
  * however in case we do want to inherit in the future in a manner
  * similar to Trajin.
  */
// FIXME: Should open state be tracked and error given if setup called
//        twice for already open traj? Should append be default if
//        traj is opened twice?
// FIXME: Should take ArgList&, not const, so we can determine if there are unknown args later.
class Trajout_Single {
  public:
    Trajout_Single() : trajio_(0), debug_(0) {}
    ~Trajout_Single();
    void SetDebug(int d) { debug_ = d; }
    // ----- Inherited functions -----------------
    /// Prepare trajectory for writing to the given format, but no Topology setup.
    int InitTrajWrite(FileName const&, ArgList const&, DataSetList const& DSLin, TrajectoryFile::TrajFormatType);
    /// Peform Topology-related setup for trajectory and open. TODO const&
    int SetupTrajWrite(Topology*, CoordinateInfo const&, int);
    /// Close output trajectory.
    void EndTraj();
    /// Write a single frame.
    int WriteSingle(int, Frame const&);
    /// Print information on trajectory to be written.
    void PrintInfo(int) const;
    // -------------------------------------------
    OutputTrajCommon Traj() const { return traj_;        }
    bool IsInitialized()    const { return trajio_ != 0; }
    /// Init and setup/open traj.
    int PrepareTrajWrite(FileName const&, ArgList const&, DataSetList const&, Topology*,
                         CoordinateInfo const&, int, TrajectoryFile::TrajFormatType);
    /// Init and setup/open traj for writing to STDOUT (e.g. ambpdb mode)
    int PrepareStdoutTrajWrite(ArgList const&, DataSetList const&, Topology*,
                               CoordinateInfo const&, int, TrajectoryFile::TrajFormatType);
    /// Init traj; if given, append ensemble number to name (use in Actions)
    int InitEnsembleTrajWrite(FileName const&, ArgList const&, DataSetList const&,
                              TrajectoryFile::TrajFormatType, int);
#   ifdef MPI
    // Set the parallel communicator.
    int SetTrajComm(Parallel::Comm const& c) { trajComm_ = c; return 0; }
#   endif
  private:
    int InitTrajout(FileName const&, ArgList const&, DataSetList const&, TrajectoryFile::TrajFormatType);
#   ifdef MPI
    /// Peform Topology-related setup for trajectory and open in parallel.
    int ParallelSetupTrajWrite();
    Parallel::Comm trajComm_;
#   endif
    OutputTrajCommon traj_;
    TrajectoryIO* trajio_;
    int debug_;
};
#endif
