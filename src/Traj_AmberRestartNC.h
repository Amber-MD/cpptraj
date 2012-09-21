#ifndef INC_TRAJ_AMBERRESTARTNC_H
#define INC_TRAJ_AMBERRESTARTNC_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
// Class: Traj_AmberRestartNC
/// TrajectoryIO class for reading and writing Amber Netcdf Restarts
class Traj_AmberRestartNC : public TrajectoryIO, NetcdfFile {
  public:
    Traj_AmberRestartNC();
    ~Traj_AmberRestartNC();
    // AmberNetcdf-specific functions
    void SetNoVelocity();

    // Inherited functions
    bool ID_TrajFormat();
    int setupTrajin(Topology*);
    int setupTrajout(Topology*,int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    int processWriteArgs(ArgList*);
    int NreplicaDimensions() { return remd_dimension_; }
    int readIndices(int,int*);
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
