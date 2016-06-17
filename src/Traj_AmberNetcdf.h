#ifndef INC_TRAJ_AMBERNETCDF_H
#define INC_TRAJ_AMBERNETCDF_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
/// Reads and writes Amber Netcdf format trajectories. 
class Traj_AmberNetcdf : public TrajectoryIO, private NetcdfFile {
  public:
    Traj_AmberNetcdf();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_AmberNetcdf(); }
    ~Traj_AmberNetcdf();
    static void ReadHelp();
    static void WriteHelp();
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&);
    int processReadArgs(ArgList&);
    // Reservoir functions
    inline int createReservoir(bool,double,int);
    int writeReservoir(int, Frame const&, double, int);
#   ifdef MPI
    // Parallel functions
    int parallelOpenTrajin(Parallel::Comm const&);
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelReadFrame(int, Frame&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj();
#   endif
  private:
    float *Coord_;
    FileName filename_;
    int eptotVID_;
    int binsVID_;
    bool useVelAsCoords_;
    bool useFrcAsCoords_;
    bool readAccess_;
    bool outputTemp_;
    bool outputVel_;
    bool outputFrc_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
int Traj_AmberNetcdf::createReservoir(bool hasBins, double reservoirT, int iseed) {
  return NC_createReservoir(hasBins, reservoirT, iseed, eptotVID_, binsVID_);
}
#endif
#endif
