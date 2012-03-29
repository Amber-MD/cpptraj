#ifndef INC_TRAJ_AMBERRESTARTNC_H
#define INC_TRAJ_AMBERRESTARTNC_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
// Class: AmberRestartNC
/// TrajectoryIO class for reading and writing Amber Netcdf Restarts
class AmberRestartNC : public TrajectoryIO {
  public:
    AmberRestartNC();
    ~AmberRestartNC();
    // AmberNetcdf-specific functions
    void SetNoVelocity();
  private:
    int ncid;
    int atomDID;
    int ncatom;
    int ncatom3;
    int coordVID;
    int velocityVID;
    double velocityScale;
    int cellAngleVID;
    int cellLengthVID;

    int spatialDID;
    int labelDID;
    int cell_spatialDID;
    int cell_angularDID;
    int spatialVID;
    int timeVID;
    double restartTime;
    int cell_spatialVID;
    int cell_angularVID;
    int TempVID;

    bool singleWrite;
    double time0;
    double dt;

    // Inherited functions
    int setupTrajin(AmberParm*);
    int setupTrajout(AmberParm*);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    int processWriteArgs(ArgList*);
    void info();

    int setupWriteForSet(int,double*);
};
#endif
#endif  
