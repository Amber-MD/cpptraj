#ifndef INC_TRAJ_AMBERRESTARTNC_H
#define INC_TRAJ_AMBERRESTARTNC_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
// Class: AmberRestartNC
/// TrajectoryIO class for reading and writing Amber Netcdf Restarts
class AmberRestartNC : public TrajectoryIO, NetcdfFile {
  public:
    AmberRestartNC();
    ~AmberRestartNC();
    // AmberNetcdf-specific functions
    void SetNoVelocity();

    // Inherited functions
    bool ID_TrajFormat();
    int setupTrajin(Topology*);
    int setupTrajout(Topology*);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    int processWriteArgs(ArgList*);
    void info();

  private:
    double restartTime_;
    bool singleWrite_;
    double time0_;
    double dt_;

    int setupWriteForSet(int,double*);
};
#endif
#endif  
