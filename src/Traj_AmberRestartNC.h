#ifndef INC_TRAJ_AMBERRESTARTNC_H
#define INC_TRAJ_AMBERRESTARTNC_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
/// Class: AmberRestartNC
/// TrajectoryIO class for reading and writing Amber Netcdf Restarts
class AmberRestartNC : public TrajectoryIO {
    int ncid;
    int atomDID;
    int ncatom;
    int ncatom3;
    int coordVID;
    bool hasVelocity;
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

    // Inherited functions
    int setupRead(int);
    int setupWrite(int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();

    int setupWriteForSet(int);

  public:
    AmberRestartNC();
    ~AmberRestartNC();
    // AmberNetcdf-specific functions
};
#endif
#endif  
