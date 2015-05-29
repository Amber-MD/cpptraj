#ifndef INC_ENSEMBLE_SINGLE_H
#define INC_ENSEMBLE_SINGLE_H
#ifdef ENABLE_SINGLE_ENSEMBLE
#include "Ensemble.h"
#include "TrajectoryIO.h"
/// Read in an array of frames at a time from a single file.
class Ensemble_Single : public Ensemble {
  public:
    Ensemble_Single();
    ~Ensemble_Single();
    // ----- Inherited Functions -----------------
    int SetupEnsembleRead(std::string const&, ArgList&, Topology*);
    int ReadEnsemble(int, FrameArray&, FramePtrArray&);
    int BeginEnsemble();
    void EndEnsemble();
    void EnsembleInfo(int) const;
    CoordinateInfo const& EnsembleCoordInfo() const { return cInfo_; }
    // -------------------------------------------
  private:
    TrajectoryIO* eio_;
    int ensembleSize_; // Should always equal cInfo_.EnsembleSize()
    CoordinateInfo cInfo_;
};
#endif
#endif
