#ifndef INC_TRAJ_AMBERNETCDF_H
#define INC_TRAJ_AMBERNETCDF_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
/// Class: AmberNetcdf
/// TrajectoryIO class for reading and writing Amber Netcdf format 
/// trajectories. 
class AmberNetcdf : public TrajectoryIO {
    int ncid;
    int frameDID;
    int ncframe;
    int atomDID;
    int ncatom; 
    int ncatom3; 
    float *Coord;
    int coordVID;
    int cellAngleVID;
    int cellLengthVID;

    int spatialDID;
    int labelDID;
    int cell_spatialDID;
    int cell_angularDID;
    int spatialVID;
    int timeVID;
    int cell_spatialVID;
    int cell_angularVID;
    int TempVID;

    // Inherited functions
    int setupRead(int);
    int setupWrite(int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*);
    int writeFrame(int,double*,double*,double);
    void info();

  public:
    AmberNetcdf();
    ~AmberNetcdf();
    // AmberNetcdf-specific functions
    void SetRemdTraj();
};
#endif
#endif
