#ifndef INC_TRAJ_AMBERNETCDF_H
#define INC_TRAJ_AMBERNETCDF_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
// Class: AmberNetcdf
/// Reads and writes Amber Netcdf format trajectories. 
class AmberNetcdf : public TrajectoryIO {
  public:
    AmberNetcdf();
    ~AmberNetcdf();
    // AmberNetcdf-specific functions
    // Inherited functions
    int setupTrajin(Topology*);
    int setupTrajout(Topology*);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList *);

  private:
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
    // Multi-D RMED
    int remd_dimension;
    int dimensionDID;
    int groupnumVID;
    int dimtypeVID;
    int indicesVID;
    int *remd_groupnum;
    int *remd_dimtype;
    int *remd_indices;

};
#endif
#endif
