#ifndef INC_TRAJ_H5MD_H
#define INC_TRAJ_H5MD_H
#include "TrajectoryIO.h"
//namespace H5 {
//  class H5File;
//}
/// MDanalysis H5MD (HDF5) format 
class Traj_H5MD : public TrajectoryIO {
  public:
    Traj_H5MD();
    ~Traj_H5MD();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_H5MD(); }
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
    /// Set up coordinates VID and related dim sizes
    int setupCoordVID(int, int&, int&, int&, int&);
    /// Set up box variable IDs and determine type
    int setupBoxVIDs(Box&, int, int);
#   endif

//    H5::H5File* file_;
    typedef std::vector<std::string> Sarray;
    typedef std::vector<int> Iarray;

    Sarray mainGroupNames_;
    Iarray mainGroupIds_;

    int particle_gid_;  ///< Particle group id
    
    int ncid_;          ///< NetCDF ID
    int natom_;         ///< Number of atoms in each trajectory frame.
    int coordVID_;      ///< Coordinates variable ID
    int cellLengthVID_; ///< Cell lengths variable ID
    int cellAngleVID_;  ///< Cell angles variable ID
    int timeVID_;       ///< Time variable ID
    size_t start_[3];   ///< Start indices in each dimension
    size_t count_[3];   ///< Count indices in each dimension
    double convert_h5_to_cpptraj_box_; ///< For converting h5 box lengths to angstroms
    double convert_h5_to_cpptraj_coord_; ///< For converting h5 coordinates to angstroms
    std::vector<float> ftmp_; ///< Temporary array to store floats
};
#endif
