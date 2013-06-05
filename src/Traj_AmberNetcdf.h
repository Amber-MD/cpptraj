#ifndef INC_TRAJ_AMBERNETCDF_H
#define INC_TRAJ_AMBERNETCDF_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
// Class: Traj_AmberNetcdf
/// Reads and writes Amber Netcdf format trajectories. 
class Traj_AmberNetcdf : public TrajectoryIO, private NetcdfFile {
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
    ReplicaDimArray ReplicaDimensions() { return remdDim_; }
    int readIndices(int,int*);
    // Reservoir functions
    inline int createReservoir(bool,double,int);
    int writeReservoir(int, Frame&, double, int);
  private:
    float *Coord_;
    FileName filename_;
    int eptotVID_;
    int binsVID_;
    int processReadArgs(ArgList&) { return 0; }
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
int Traj_AmberNetcdf::createReservoir(bool hasBins, double reservoirT, int iseed) {
  return NC_createReservoir(hasBins, reservoirT, iseed, eptotVID_, binsVID_);
}
#endif
#endif
