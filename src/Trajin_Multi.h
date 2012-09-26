#ifndef INC_TRAJIN_MULTI_H
#define INC_TRAJIN_MULTI_H
#include <map>
#include "Trajin.h"
#include "FrameArray.h"
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

    int EnsembleSetup( FrameArray& );
    int GetNextEnsemble( FrameArray& );
    int EnsembleSize() { return (int)REMDtraj_.size(); }
    int EnsemblePosition(int member) { return frameidx_[member]; }
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
    TargetType targetType_;   ///< Hold type of REMD frame being searched for.
    NameListType replica_filenames_;
    // ENSEMBLE
    RemdIdxType frameidx_;    ///< Hold position of each frame in ensemble.
    typedef std::map<double,int> TmapType;
    TmapType TemperatureMap_;

    NameListType SearchForReplicas(bool);
    bool IsTarget(double);
};
#endif
