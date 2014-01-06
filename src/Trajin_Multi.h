#ifndef INC_TRAJIN_MULTI_H
#define INC_TRAJIN_MULTI_H
#include <map>
#include "Trajin.h"
#include "FrameArray.h"
#include "DataSet_RemLog.h"
/// Class for reading in multiple trajectories at the same time (e.g. REMD ensemble)
class Trajin_Multi : public Trajin {
  public:
    Trajin_Multi();
    ~Trajin_Multi();

    int SetupTrajRead(std::string const&, ArgList&, Topology*);
    int BeginTraj(bool);
    void EndTraj();
    int ReadTrajFrame( int, Frame& );
    void PrintInfo(int) const;
    bool HasVelocity()      const { return hasVelocity_; }
    int NreplicaDimension() const { return Ndimensions_; }

    void EnsembleInfo() const;
    int EnsembleSetup( FrameArray& );
    int GetNextEnsemble( FrameArray& );
    int EnsembleSize()               const { return (int)REMDtraj_.size(); }
#   ifdef MPI
    int EnsembleFrameNum()           const { return ensembleFrameNum_;     }
#   else
    int EnsemblePosition(int member) const { return frameidx_[member];     }
#   endif
    bool BadEnsemble()               const { return badEnsemble_;          }
    // CRDIDXARG: NOTE: This is public for CRDIDX in TrajinList
    enum TargetType { NONE = 0, TEMP, INDICES, CRDIDX };
    TargetType TargetMode()          const { return targetType_;           }
    std::string FinalCrdIndices()    const;
  private:
    /// Define type that will hold REMD indices
    typedef std::vector<int> RemdIdxType;
    typedef std::vector<TrajectoryIO*> IOarrayType;
    typedef std::vector<std::string> NameListType;

    double remdtrajtemp_;     ///< Get frames with this temperature on read
    RemdIdxType remdtrajidx_; ///< Get frames with these indices on read
    int Ndimensions_;         ///< # of dimensions in each trajectory.
    IOarrayType REMDtraj_;    ///< Input replica trajectories
    int lowestRepnum_;        ///< Hold the lowest replica number
    bool hasVelocity_;        ///< True if all trajs have velocities.
    bool replicasAreOpen_;    ///< True is replicas are open.
    bool badEnsemble_;        ///< True if problem with any frames in the ensemble
    TargetType targetType_;   ///< Hold type of REMD frame being searched for.
    NameListType replica_filenames_;
    // ENSEMBLE
    //RemdIdxType frameidx_;    ///< Hold position of each frame in ensemble.
    int* frameidx_;    ///< Hold position of each frame in ensemble.
    typedef std::map<double,int> TmapType;
    TmapType TemperatureMap_;
    typedef std::map< std::vector<int>, int > ImapType;
    ImapType IndicesMap_;
#   ifdef MPI
    int ensembleFrameNum_;      ///< Position containing coords to use in FrameArray
#   endif
    DataSet_RemLog remlogData_; ///< For sorting by CRDIDX from remlog.
    NameListType SearchForReplicas();
    bool IsTarget(Frame const&);
};
#endif
