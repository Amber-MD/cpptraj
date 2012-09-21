#ifndef INC_REMDTRAJ_H
#define INC_REMDTRAJ_H
#include "TrajectoryIO.h"
#include <vector>
// Class: RemdTraj
/// Used to process input trajectories in parallel.
/** The primary use is to read in replica trajectories. Hold each replica 
  * trajectory in its own TrajectoryIO. During reads, only the frame with 
  * a temperature matching remdtrajtemp will be used.
  */
class RemdTraj : public TrajectoryIO {
  public:
    typedef std::vector<int> RemdIdxType;

    RemdTraj();
    ~RemdTraj();
    
    // RemdTraj-specific functions
    void SetTargetTemp(double);
    void SetTargetIdx(RemdIdxType const&);
    int SetupTemperatureList(int);
    std::vector<std::string> SearchForReplicas();
    int GetTemperatureName(std::string&,const char*, int);
    void AddReplicaTrajin(TrajectoryIO*);
    void AddReplicaTrajout(TrajectoryIO*);
    const char *LowestReplicaName();

    int Nreplicas();
  private:
    enum TargetType { TEMP = 0, INDICES };
    // RemdTraj-specific variables 
    double remdtrajtemp_;                    ///< Get frames with this temperature on read
    RemdIdxType remdtrajidx_;                ///< Get frames with these indices on read
    int* remd_indices_;                      ///< Space for reading in REMD indices.
    std::vector<TrajectoryIO*> REMDtraj_;    ///< Input replica trajectories
    std::vector<TrajectoryIO*> REMDtrajout_; ///< Output replica trajectories
    int lowestRepnum_;                       ///< Hold the lowest replica number
    TargetType targetType_;

    // For T-trajectory writes 
    bool hasTrajout_;         ///< True if writing replica trajectories during read
    double *remdX_;           ///< Hold coords of replica traj for writing out
    double *remdV_;           ///< Hold velocities of rep traj for writing out
    double remdbox_[6];       ///< Hold box coords of rep traj for writing out
    double remdT_;            ///< Hold replica temp for writing out 
    int remdN_;               ///< For output trajectories, number of coords.
    std::vector<double> TemperatureList_; ///< List of temperatures found in replicas

    void PrintNoExtError();
    bool IsTarget(double);

    // Inherited functions
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    //int writeFrame(int,double*,double*,double);
    void info();
};
#endif
