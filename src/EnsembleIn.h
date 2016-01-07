#ifndef INC_ENSEMBLEIN_H
#define INC_ENSEMBLEIN_H
#include "InputTrajCommon.h"
#include "FrameArray.h"
#include "FramePtrArray.h"
#include "ReplicaInfo.h"
#ifdef MPI
# include "Parallel.h"
# ifdef TIMER
#   include "Timer.h"
# endif
#endif
/// Read in an array of frames at a time.
class EnsembleIn {
  public:
    EnsembleIn() : targetType_(ReplicaInfo::NONE), badEnsemble_(0), debug_(0) {}
    virtual ~EnsembleIn() {}
    virtual int SetupEnsembleRead(FileName const&, ArgList&, Topology*) = 0;
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
    /// \return Common information (counter, filename, Topology ptr)
    InputTrajCommon const& Traj() const { return traj_;        }
    /// \return true if there was a problem reading the ensemble.
    int BadEnsemble()             const { return badEnsemble_; }
    void SetDebug(int d)                { debug_ = d;          }
    /// Write current replica mapping to STDOUT.
    void PrintReplicaInfo() const;
  protected:
    typedef Frame::RemdIdxType RemdIdxType;

    int SetTemperatureMap(std::vector<double> const&);
    int SetIndicesMap(std::vector<RemdIdxType> const&);
    InputTrajCommon& SetTraj() { return traj_; }
    /// For converting temperature to replica index
    ReplicaMap<double> TemperatureMap_;
    /// For converting indices to replica index
    ReplicaMap<RemdIdxType> IndicesMap_;
    ReplicaInfo::TargetType targetType_; ///< Hold type of REMD frame being searched for.
    int badEnsemble_;         ///< Set to 1 if problem reading ensemble, 0 otherwise.
    int debug_;
#   ifdef MPI
    RemdIdxType frameidx_;    ///< Hold position of each frame in ensemble.

    int Member()                         const { return Parallel::EnsembleComm().Rank(); }
    Parallel::Comm const& EnsembleComm() const { return Parallel::EnsembleComm();        }
    static int GatherTemperatures(double*, std::vector<double>&, Parallel::Comm const&);
    static int GatherIndices(int*, std::vector<RemdIdxType>&, int, Parallel::Comm const&);
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
int EnsembleIn::GetNextEnsemble(FrameArray& fa, FramePtrArray& fp) {
  if (traj_.Counter().CheckFinished()) return 0;
  if (ReadEnsemble( traj_.Counter().Current(), fa, fp )) return 0;
  traj_.Counter().UpdateCounters();
  return 1;
}
#endif
