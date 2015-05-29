#ifndef INC_ENSEMBLE_H
#define INC_ENSEMBLE_H
#include "InputTrajCommon.h"
#include "FrameArray.h"
#include "FramePtrArray.h"
#include "ReplicaInfo.h"
#ifdef MPI
#  ifdef TIMER
#    include "Timer.h"
#  endif
#endif
/// Read in an array of frames at a time.
class Ensemble {
  public:
    Ensemble() : targetType_(ReplicaInfo::NONE), badEnsemble_(0), debug_(0) {}
    virtual ~Ensemble() {}
    virtual int SetupEnsembleRead(std::string const&, ArgList&, Topology*) = 0;
    virtual int ReadEnsemble(int, FrameArray&, FramePtrArray&) = 0;
    virtual int BeginEnsemble() = 0;
    virtual void EndEnsemble() = 0;
    virtual void EnsembleInfo(int) const = 0;
    virtual CoordinateInfo const& EnsembleCoordInfo() const = 0;
    // -------------------------------------------
#   ifdef MPI
#   ifdef TIMER
    static void TimingData(double);
#   endif
#   endif
    inline int GetNextEnsemble(FrameArray&, FramePtrArray&);

    InputTrajCommon const& Traj() const { return traj_;        }
    /// \return true if there was a problem reading the ensemble.
    int BadEnsemble()             const { return badEnsemble_; }
    void SetDebug(int d)                { debug_ = d;          }
    /// Write current replica mapping to STDOUT.
    void PrintReplicaInfo() const;
  protected:
    typedef Frame::RemdIdxType RemdIdxType;
    InputTrajCommon& Traj()             { return traj_; }

    int SetTemperatureMap(std::vector<double> const&);
    int SetIndicesMap(std::vector<RemdIdxType> const&);

    /// For converting temperature to replica index
    ReplicaMap<double> TemperatureMap_;
    /// For converting indices to replica index
    ReplicaMap<RemdIdxType> IndicesMap_;
    ReplicaInfo::TargetType targetType_; ///< Hold type of REMD frame being searched for.
    int badEnsemble_;         ///< Set to 1 if problem reading ensemble, 0 otherwise.
    int debug_;
#   ifdef MPI
    RemdIdxType frameidx_;    ///< Hold position of each frame in ensemble.

    int GatherTemperatures(const double*, std::vector<double>&);
    int GatherIndices(const double*, std::vector<int>&, int);
#   ifdef TIMER
    Timer mpi_allgather_timer_;
    Timer mpi_sendrecv_timer_;
    static double total_mpi_allgather_;
    static double total_mpi_sendrecv_;
#   endif
#   endif
  private:
    static void PrintReplicaTmap(ReplicaMap<double> const&);
    static void PrintReplicaImap(ReplicaMap<RemdIdxType> const&);

    InputTrajCommon traj_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
int Ensemble::GetNextEnsemble(FrameArray& fa, FramePtrArray& fp) {
  if (traj_.Counter().CheckFinished()) return 0;
  if (ReadEnsemble( traj_.Counter().Current(), fa, fp )) return 0;
  traj_.Counter().UpdateCounters();
  return 1;
}
#endif
