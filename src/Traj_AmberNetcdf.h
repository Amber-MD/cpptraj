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
    int processWriteArgs(ArgList&, DataSetList const&);
    int processReadArgs(ArgList&);
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
    float *Coord_;        ///< Temporary array for converting double <-> single precision
    FileName filename_;   ///< File name
    bool useVelAsCoords_; ///< If true read velocities in place of coordinates
    bool useFrcAsCoords_; ///< If true read forces in place of coordinates
    bool readAccess_;     ///< True if set up for read, false otherwise
    bool outputTemp_;     ///< If true write out temperatures
    bool write_mdcrd_;    ///< If true write out coordinates
    bool write_mdvel_;    ///< If true write out velocities
    bool write_mdfrc_;    ///< If true write out forces
};
#endif
#endif
