#ifndef INC_ENSEMBLEOUT_MULTI_H
#define INC_ENSEMBLEOUT_MULTI_H
#include "EnsembleOut.h"
/// Write out an array of frames at a time to multiple files.
class EnsembleOut_Multi : public EnsembleOut {
  public:
    EnsembleOut_Multi();
    ~EnsembleOut_Multi();
    // ----- Inherited Functions -----------------
    int InitEnsembleWrite(std::string const&, ArgList const&, DataSetList const&, int, TrajectoryFile::TrajFormatType);
    int SetupEnsembleWrite(Topology*, CoordinateInfo const&, int);
    void EndEnsemble();
    int WriteEnsemble(int, FramePtrArray const&);
    void PrintInfo(int) const;
  private:
    void Clear();
#   ifdef MPI
    int ParallelSetupEnsembleWrite();
#   endif
    typedef std::vector<TrajectoryIO*> IOarrayType;
    IOarrayType ioarray_;
    typedef std::vector<std::string> Sarray;
    Sarray fileNames_;
    int ensembleSize_;
#   ifndef MPI
    std::vector<int> tIndex_;
#   endif
};
#endif
