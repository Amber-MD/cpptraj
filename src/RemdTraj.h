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
    RemdTraj();
    ~RemdTraj();
    
    // RemdTraj-specific functions
    void SetTargetTemp(double);
    int SetupTemperatureList(int);
    std::vector<std::string> SearchForReplicas();
    int GetTemperatureName(std::string&,char*, int);
    void AddReplicaTrajin(TrajectoryIO*);
    void AddReplicaTrajout(TrajectoryIO*);
    const char *LowestReplicaName();

    int Nreplicas();
  private:
    // RemdTraj-specific variables 
    double remdtrajtemp_;                    ///< Get frames with this temperature on read
    std::vector<TrajectoryIO*> REMDtraj_;    ///< Input replica trajectories
    std::vector<TrajectoryIO*> REMDtrajout_; ///< Output replica trajectories
    int lowestRepnum_;                       ///< Hold the lowest replica number

    // For T-trajectory writes 
    bool hasTrajout_;         ///< True if writing replica trajectories during read
    double *remdX_;           ///< Hold coords of replica traj for writing out
    double *remdV_;           ///< Hold velocities of rep traj for writing out
    double remdbox_[6];       ///< Hold box coords of rep traj for writing out
    double remdT_;            ///< Hold replica temp for writing out 
    int remdN_;               ///< For output trajectories, number of coords.
    std::vector<double> TemperatureList_; ///< List of temperatures found in replicas

    void PrintNoExtError();

    // Inherited functions
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    //int writeFrame(int,double*,double*,double);
    void info();
};
#endif
