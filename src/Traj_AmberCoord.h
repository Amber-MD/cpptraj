#ifndef INC_TRAJ_AMBERCOORD_H
#define INC_TRAJ_AMBERCOORD_H
#include "TrajectoryIO.h"
#include "BufferedFile.h"
// Class: Traj_AmberCoord
/// Reads and writes formatted (ASCII text) amber trajectories. 
class Traj_AmberCoord: public TrajectoryIO {
  public:
    Traj_AmberCoord();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_AmberCoord(); }
    ~Traj_AmberCoord();
  private:
    int debug_;           ///< Debug level
    int natom3_;          ///< Number of coords (# atoms X 3)
    size_t titleSize_;    ///< Title size in bytes
    size_t hasREMD_;      ///< REMD header size if present
    int numBoxCoords_;    ///< Number of box coords, 3 (ortho or truncoct) or 6 (triclinic)
    const char* outfmt_;  ///< Format string for writing coordinates
    bool highPrecision_;  ///< If true output format will be 8.6 instead of 8.3
    TrajectoryIO* mdvel_; ///< Associated file with velocities.
    BufferedFile file_;   ///< Buffer reads/writes
    std::string title_;   ///< Trajectory title (max 80 chars)
    bool seekable_;       ///< True if traj can be randomly accessed.
    double boxAngle_[3];  ///< Hold default box angles in case traj has only box lengths

    static const size_t REMD_HEADER_SIZE;
    static const size_t BUF_SIZE;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*, TrajInfo&);
    int setupTrajout(std::string const&, Topology*, int, TrajInfo const&, bool);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList&);
    int processReadArgs(ArgList&);
    int readVelocity(int, double*) { return 1; }
    int readIndices(int,int*) { return 1; }
};
#endif
