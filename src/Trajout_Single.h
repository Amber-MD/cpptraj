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
    int WriteSingle(int, Frame const&);
    int WriteEnsemble(int,FramePtrArray const&) { return 1; }
    void PrintInfo(int) const;
    int SetupTrajWrite(Topology*);
    // -------------------------------------------
    /// For writing single traj to STDOUT (e.g. ambpdb mode)
    int PrepareStdoutTrajWrite(ArgList const&, Topology*, TrajectoryFile::TrajFormatType);
    /// For writing single traj from Action, ensemble-aware.
    int InitEnsembleTrajWrite(std::string const&, ArgList const&, Topology*,
                              TrajFormatType fmtIn, int ensembleNum);
  private:
    int InitTrajout(std::string const&, ArgList const&, Topology*, TrajectoryFile::TrajFormatType);

    TrajectoryIO* trajio_;
};
#endif
