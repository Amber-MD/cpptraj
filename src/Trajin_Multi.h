#ifndef INC_TRAJIN_MULTI_H
#define INC_TRAJIN_MULTI_H
#include "Trajin.h"
#include "DataSet_RemLog.h"
#ifdef MPI
#  ifdef TIMER
#    include "Timer.h"
#  endif
#endif
/// Class for reading in multiple trajectories at the same time (e.g. REMD ensemble)
class Trajin_Multi : public Trajin {
  public:
    Trajin_Multi();
    ~Trajin_Multi();
    // ----- Inherited functions -----------------
    int SetupTrajRead(std::string const&, ArgList&, Topology*);
    int ReadTrajFrame( int, Frame& );
    int BeginTraj(bool);
    void EndTraj();
    void PrintInfo(int) const;
    CoordinateInfo const& TrajCoordInfo() const { return cInfo_; }

    // NOTE: The following are currently for testing Trajin_Ensemble
    void EnsembleInfo() const;
    int EnsembleSetup( FrameArray&, FramePtrArray& );
    /// \return 0 if more frames to read, 1 if finished/error
    int ReadEnsemble( int, FrameArray&, FramePtrArray& );
    bool BadEnsemble()      const { return badEnsemble_;          }
    // -------------------------------------------
    // CRDIDXARG
    ReplicaInfo::TargetType TargetMode() const { return targetType_; }
    std::string FinalCrdIndices()        const; // TODO: Obsolete.
#   ifdef MPI
#   ifdef TIMER
    static void TimingData(double);
#   endif
#   endif
  private:
    /// Define type that will hold REMD indices
    typedef Frame::RemdIdxType RemdIdxType;
    typedef std::vector<TrajectoryIO*> IOarrayType;
    typedef std::vector<std::string> NameListType;

    NameListType SearchForReplicas();
    bool IsTarget(Frame const&);

    CoordinateInfo cInfo_;    ///< Collective coord information for all replicas.
    double remdtrajtemp_;     ///< Get frames with this temperature on read
    double remdFrameFactor_;  ///< For HREMD sort, # frames written per remlog entry
    int remdFrameOffset_;     ///< If traj written less often than log, +1
    RemdIdxType remdtrajidx_; ///< Get frames with these indices on read
    IOarrayType REMDtraj_;    ///< Input replica trajectories
    int lowestRepnum_;        ///< Hold the lowest replica number
    bool replicasAreOpen_;    ///< True is replicas are open.
    bool badEnsemble_;        ///< True if problem with any frames in the ensemble
    ReplicaInfo::TargetType targetType_; ///< Hold type of REMD frame being searched for.
    NameListType replica_filenames_;
    // ENSEMBLE
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
    DataSet_RemLog remlogData_; ///< For sorting by CRDIDX from remlog.
};
#endif
