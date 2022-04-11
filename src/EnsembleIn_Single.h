#ifndef INC_ENSEMBLEIN_SINGLE_H
#define INC_ENSEMBLEIN_SINGLE_H
#ifdef ENABLE_SINGLE_ENSEMBLE
#include "EnsembleIn.h"
// Forward declarations
class TrajectoryIO;
/// Read in an array of frames at a time from a single file.
class EnsembleIn_Single : public EnsembleIn {
  public:
    EnsembleIn_Single();
    ~EnsembleIn_Single();
    // ----- Inherited Functions -----------------
    int SetupEnsembleRead(FileName const&, ArgList&, Topology*);
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
