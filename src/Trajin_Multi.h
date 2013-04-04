#ifndef INC_TRAJIN_MULTI_H
#define INC_TRAJIN_MULTI_H
#include <map>
#include "Trajin.h"
#include "FrameArray.h"
/// Class for reading in multiple trajectories at the same time (e.g. REMD ensemble)
class Trajin_Multi : public Trajin {
  public:
    Trajin_Multi();
    ~Trajin_Multi();

    int SetupTrajRead(std::string const&, ArgList *, Topology *);
    int BeginTraj(bool);
    void EndTraj();
    int GetNextFrame( Frame& );
    void PrintInfo(int);
    bool HasVelocity() { return hasVelocity_; }

    void EnsembleInfo() const;
    int EnsembleSetup( FrameArray& );
    int GetNextEnsemble( FrameArray& );
    int EnsembleSize()               const { return (int)REMDtraj_.size(); }
    int EnsemblePosition(int member) const { return frameidx_[member];     }
    bool BadEnsemble()               const { return badEnsemble_;          }
  private:
    /// Define type that will hold REMD indices
    typedef std::vector<int> RemdIdxType;
    typedef std::vector<TrajectoryIO*> IOarrayType;
    typedef std::vector<std::string> NameListType;
    enum TargetType { NONE = 0, TEMP, INDICES };

    double remdtrajtemp_;     ///< Get frames with this temperature on read
    RemdIdxType remdtrajidx_; ///< Get frames with these indices on read
    int* remd_indices_;       ///< Space for reading in REMD indices.
    IOarrayType REMDtraj_;    ///< Input replica trajectories
    int lowestRepnum_;        ///< Hold the lowest replica number
    bool isSeekable_;         ///< True if all trajs are seekable.
    bool hasVelocity_;        ///< True if all trajs have velocities.
    bool isEnsemble_;         ///< True if this will be processed as an ensemble.
    bool replicasAreOpen_;    ///< True is replicas are open.
    bool badEnsemble_;        ///< True if problem with any frames in the ensemble
    TargetType targetType_;   ///< Hold type of REMD frame being searched for.
    NameListType replica_filenames_;
    // ENSEMBLE
    RemdIdxType frameidx_;    ///< Hold position of each frame in ensemble.
    typedef std::map<double,int> TmapType;
    TmapType TemperatureMap_;

    NameListType SearchForReplicas();
    bool IsTarget(double);
};
#endif
