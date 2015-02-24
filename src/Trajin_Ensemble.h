#ifndef TRAJIN_ENSEMBLE_H
#define TRAJIN_ENSEMBLE_H
#ifdef ENABLE_SINGLE_ENSEMBLE
#include "Trajin.h"
#ifdef MPI
#  ifdef TIMER
#    include "Timer.h"
#  endif
#endif
/// Class for reading in single file ensemble trajectories.
class Trajin_Ensemble : public Trajin {
  public:
    Trajin_Ensemble();
    ~Trajin_Ensemble();
    // ----- Inherited Functions -----------------
    int SetupTrajRead(std::string const&, ArgList&, Topology*);
    int ReadTrajFrame(int, Frame&) {return 1;}
    int BeginTraj(bool);
    void EndTraj();
    void PrintInfo(int) const;
    CoordinateInfo const& TrajCoordInfo() const { return cInfo_; }

    // NOTE: The following are currently for testing Trajin_Ensemble
    void EnsembleInfo() const;
    int EnsembleSetup(FrameArray&, FramePtrArray&);
    int ReadEnsemble(int, FrameArray&, FramePtrArray&);
    bool BadEnsemble() const { return badEnsemble_; }
    // -------------------------------------------
#   ifdef MPI
#   ifdef TIMER
    static void TimingData(double);
#   endif
#   endif
  private:
    typedef Frame::RemdIdxType RemdIdxType;
    ReplicaInfo::TargetType targetType_;
    TrajectoryIO* eio_;
    int ensembleSize_; // Should always equal cInfo_.EnsembleSize()
    bool trajIsOpen_;
    bool badEnsemble_;
    CoordinateInfo cInfo_; 
    ReplicaMap<double> TemperatureMap_;
    ReplicaMap<RemdIdxType> IndicesMap_;
#   ifdef MPI
    RemdIdxType frameidx_;    ///< Hold position of each frame in ensemble.
#   ifdef TIMER
    Timer mpi_allgather_timer_;
    Timer mpi_sendrecv_timer_;
    static double total_mpi_allgather_;
    static double total_mpi_sendrecv_;
#   endif
#   endif 
};
#endif
#endif
