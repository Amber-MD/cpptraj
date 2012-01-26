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
    // Private vars
    char *Prefix;      ///< Complete filename to lowest replica up to the numerical extension
    int ExtWidth;      ///< Size of the numerical extension in characters
    char *CompressExt; ///< If replica is compressed, hold the compression extension
    char *repFilename; ///< Will hold replica filename last set by GetReplicaName
    int lowestRepnum;  ///< Hold the lowest replica number
    bool hasTrajout;   ///< True if writing replica trajectories during read
    double *remdX;     ///< Hold coords of replica traj for writing out
    double *remdV; 
    double remdbox[6];
    double remdT;
    int remdN;         ///< For output trajectories, number of coords.
    char *repOutname;  ///< Will hold output temperature trajectory filename
    double *TemperatureList; ///< List of temperatures found in replicas

    // Inherited functions
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    //int writeFrame(int,double*,double*,double);
    void info();

  public:
    RemdTraj();
    ~RemdTraj();
    
    // RemdTraj-specific variables 
    double remdtrajtemp;                    ///< Get frames with this temperature on read
    std::vector<TrajectoryIO*> REMDtraj;    ///< Input replica trajectories
    std::vector<TrajectoryIO*> REMDtrajout; ///< Output replica trajectories
    
    // RemdTraj-specific functions
    int SetupTemperatureList(int);
    int SetReplicaName(char*);
    char *GetReplicaName(int);
    char *GetLowestReplicaName();
    char *GetTemperatureName(char *, int);
};
#endif
