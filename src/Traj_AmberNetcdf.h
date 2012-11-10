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
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_AmberNetcdf(); }
    ~Traj_AmberNetcdf();
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int readVelocity(int, double*);
    int writeFrame(int,double*,double*,double*,double);
    void Info();
    int processWriteArgs(ArgList&);
    int NreplicaDimensions() { return remd_dimension_; }
    int readIndices(int,int*);

  private:
    float *Coord_;
    float *Veloc_;
    FileName filename_;
    int processReadArgs(ArgList&) { return 0; }
};
#endif
#endif
