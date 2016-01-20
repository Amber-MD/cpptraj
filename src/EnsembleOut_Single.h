#ifndef INC_ENSEMBLEOUT_SINGLE_H
#define INC_ENSEMBLEOUT_SINGLE_H
#ifdef ENABLE_SINGLE_ENSEMBLE
#include "EnsembleOut.h"
/// Write out and array of frames at a time to single file.
class EnsembleOut_Single : public EnsembleOut {
  public:
    EnsembleOut_Single();
    ~EnsembleOut_Single();
    // ----- Inherited Functions -----------------
    int InitEnsembleWrite(std::string const&, ArgList const&, int, TrajectoryFile::TrajFormatType);
    int SetupEnsembleWrite(Topology*, CoordinateInfo const&, int);
    void EndEnsemble();
    int WriteEnsemble(int, FramePtrArray const&);
    void PrintInfo(int) const;
#   ifdef MPI
    int ParallelSetupEnsembleWrite();
#   endif
  private:
    TrajectoryIO* eio_; // TODO Make EnsembleIO
    int ensembleSize_;
};
#endif
#endif
