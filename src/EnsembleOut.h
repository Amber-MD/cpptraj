#ifndef INC_ENSEMBLEOUT_H
#define INC_ENSEMBLEOUT_H
#include "OutputTrajCommon.h"
/// Write out an array of frames at a time.
class EnsembleOut {
  public:
    EnsembleOut() : debug_(0) {}
    virtual ~EnsembleOut() {} // inherited
    void SetDebug(int d) { debug_ = d; }
    // ----- Inherited Functions -----------------
    /// Prepare ensemble for writing given format, no Topology setup.
    virtual int InitEnsembleWrite(std::string const&, ArgList const&,
                                  int, TrajectoryFile::TrajFormatType) = 0;
    /// Perform Topology-related setup for ensemble and open. TODO const&
    virtual int SetupEnsembleWrite(Topology*, CoordinateInfo const&, int) = 0;
    /// Close output ensemble.
    virtual void EndEnsemble() = 0;
    /// Write array of frames.
    virtual int WriteEnsemble(int, FramePtrArray const&) = 0;
    /// Print information on ensemble to be written
    virtual void PrintInfo(int) const = 0;
    // -------------------------------------------
    OutputTrajCommon const& Traj() const { return traj_; }
#   ifdef MPI
    // Set the parallel communicator.
    int SetTrajComm(Parallel::Comm const& c) { trajComm_ = c; return 0; }
#   endif
  protected:
    OutputTrajCommon& SetTraj() { return traj_; }
    /// For ensemble trajouts, get range of members to write.
    int SetMembersToWrite(std::string const&,int);
    Range const& MembersToWrite() const { return members_to_write_; }
    int debug_;
#   ifdef MPI
    /// Peform Topology-related setup for ensemble and open in parallel.
    virtual int ParallelSetupEnsembleWrite() = 0;
    Parallel::Comm trajComm_; //TODO: make private
#   endif
  private:
    OutputTrajCommon traj_;
    Range members_to_write_; ///< Range of members to write to.
};
#endif
