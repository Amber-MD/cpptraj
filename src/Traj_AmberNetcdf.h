#ifndef INC_TRAJ_AMBERNETCDF_H
#define INC_TRAJ_AMBERNETCDF_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
// Class: Traj_AmberNetcdf
/// Reads and writes Amber Netcdf format trajectories. 
class Traj_AmberNetcdf : public TrajectoryIO, NetcdfFile {
  public:
    Traj_AmberNetcdf();
    ~Traj_AmberNetcdf();
    // AmberNetcdf-specific functions
    // Inherited functions
    bool ID_TrajFormat();
    int setupTrajin(Topology*);
    int setupTrajout(Topology*,int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList *);

  private:
    float *Coord_;
    float *Veloc_;

    // Multi-D RMED
    int remd_dimension_;
    int dimensionDID_;
    int groupnumVID_;
    int dimtypeVID_;
    int indicesVID_;
    int *remd_groupnum_;
    int *remd_dimtype_;
    int *remd_indices_;
};
#endif
#endif
