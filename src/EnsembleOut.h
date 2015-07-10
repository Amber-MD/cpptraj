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
    virtual int SetupEnsembleWrite(Topology*) = 0;
    /// Close output ensemble.
    virtual void EndEnsemble() = 0;
    /// Write array of frames.
    virtual int WriteEnsemble(int, FramePtrArray const&) = 0;
    /// Print information on ensemble to be written
    virtual void PrintInfo(int) const = 0;
  protected:
    /// For ensemble trajouts, get range of members to write.
    static Range MembersToWrite(std::string const&,int);

    int debug_;
};
#endif
