#ifndef INC_TRAJ_H5_H
#define INC_TRAJ_H5_H
#include "TrajectoryIO.h"
//namespace H5 {
//  class H5File;
//}
/// MDtraj H5 (HDF5) format 
class Traj_H5 : public TrajectoryIO {
  public:
    Traj_H5();
    ~Traj_H5();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_H5(); }
    static void WriteHelp();
    static void ReadHelp();
  private:
    // ----- Inherited functions -----------------
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int readVelocity(int, Frame&);
    int readForce(int, Frame&);
    int processWriteArgs(ArgList&, DataSetList const&);
    int processReadArgs(ArgList&);
    // -------------------------------------------
#   ifdef MPI
    // ----- Parallel functions ------------------
    int parallelOpenTrajin(Parallel::Comm const&);
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelReadFrame(int, Frame&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj();
    // -------------------------------------------
#   endif
#   ifdef HAS_HDF5
    static bool HasConventions(int);
    /// Set up coordinates VID and related dim sizes
    int setupCoordVID(int&, int&, int&, int&);
#   endif

//    H5::H5File* file_;
    int ncid_;     ///< NetCDF ID
    int natom_;    ///< Number of atoms in each trajectory frame.
    int coordVID_; ///< Coordinates variable ID
};
#endif
