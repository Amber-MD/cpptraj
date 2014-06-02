#ifndef INC_TRAJOUT_SINGLE_H
#define INC_TRAJOUT_SINGLE_H
#include "Trajout.h"
/// Single file output trajectory class.
class Trajout_Single : public Trajout {
  public:
    Trajout_Single();
    ~Trajout_Single();
    // ----- Inherited functions -----------------
    int InitTrajWrite(std::string const&, ArgList const&, Topology*, 
                      TrajectoryFile::TrajFormatType);
    void EndTraj();
    int WriteFrame(int, Topology*, Frame const&);
    int WriteEnsemble(int,Topology*,FrameArray const&,Frame::RemdIdxType const&) { return 1; }
    void PrintInfo(int) const;
    void SetEnsembleInfo(int) {}
    // -------------------------------------------
    int InitStdoutTrajWrite(ArgList const&, Topology*, TrajectoryFile::TrajFormatType);
  private:
    int InitTrajout(std::string const&, ArgList const&, Topology*, TrajectoryFile::TrajFormatType);

    TrajectoryIO* trajio_;
};
#endif
