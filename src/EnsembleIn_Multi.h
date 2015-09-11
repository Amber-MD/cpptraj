#ifndef INC_ENSEMBLEIN_MULTI_H
#define INC_ENSEMBLEIN_MULTI_H
#include "EnsembleIn.h"
#include "TrajIOarray.h"
#include "DataSet_RemLog.h" // TODO remove
/// Read in an array of frames at a time from multiple files.
class EnsembleIn_Multi : public EnsembleIn {
  public:
    EnsembleIn_Multi() : remdFrameFactor_(1.0), remdFrameOffset_(0) {}
    ~EnsembleIn_Multi() { EndEnsemble(); }
    // ----- Inherited Functions -----------------
    int SetupEnsembleRead(FileName const&, ArgList&, Topology*);
    int ReadEnsemble(int, FrameArray&, FramePtrArray&);
    int BeginEnsemble();
    void EndEnsemble();
    void EnsembleInfo(int) const;
    CoordinateInfo const& EnsembleCoordInfo() const { return cInfo_; }
    // -------------------------------------------
    // CRDIDXARG
    std::string FinalCrdIndices() const; // TODO Remove
    ReplicaInfo::TargetType TargetMode() const { return targetType_; } // TODO Remove
  private:
    TrajIOarray REMDtraj_;
    CoordinateInfo cInfo_;
    DataSet_RemLog remlogData_; ///< For sorting by CRDIDX from remlog.
    double remdFrameFactor_;  ///< For HREMD sort, # frames written per remlog entry
    int remdFrameOffset_;     ///< If traj written less often than log, +1
};
#endif
