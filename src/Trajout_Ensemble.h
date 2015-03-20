#ifndef INC_TRAJOUT_ENSEMBLE_H
#define INC_TRAJOUT_ENSEMBLE_H
#ifdef ENABLE_SINGLE_ENSEMBLE
#include "Trajout.h"
/// Class for writing ensemble as single file.
class Trajout_Ensemble : public Trajout {
  public:
    Trajout_Ensemble();
    ~Trajout_Ensemble();
    // ----- Inherited functions -----------------
    inline int InitTrajWrite(std::string const&, ArgList const&, Topology*,
                             TrajectoryFile::TrajFormatType);
    void EndTraj();
    int WriteSingle(int, Frame const&) { return 1; }
    int WriteEnsemble(int,FramePtrArray const&);
    int SetupTrajWrite(Topology*);
    void PrintInfo(int) const;
    // -------------------------------------------
  private:
    TrajectoryIO* eio_;
    int ensembleSize_;
};
#endif
#endif
