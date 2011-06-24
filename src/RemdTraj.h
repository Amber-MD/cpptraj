#ifndef INC_REMDTRAJ_H
#define INC_REMDTRAJ_H
/// Class: RemdTraj
/// Used to read in replica trajectories. Hold each replica trajectory
/// in its own TrajectoryIO. During reads, only the frame with a temperature
/// matching remdtrajtemp will be used.
#include <vector>
class RemdTraj : public TrajectoryIO {

  public:
    double remdtrajtemp;
    std::vector<TrajectoryIO*> REMDtraj;
    std::vector<TrajectoryIO*> REMDtrajout;
    double *TemperatureList;
    int numReplicas;


    RemdTraj();
    ~RemdTraj();
    // Inherited functions
    int openTraj();
    void closeTraj();
    int readFrame(int,double*, double*,double*);
    void info(int);
};
#endif
