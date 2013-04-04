#ifndef INC_TRAJ_AMBERCOORD_H
#define INC_TRAJ_AMBERCOORD_H
#include "TrajectoryIO.h"
#include "BufferedFrame.h"
// Class: Traj_AmberCoord
/// Reads and writes formatted (ASCII text) amber trajectories. 
class Traj_AmberCoord: public TrajectoryIO {
  public:
    Traj_AmberCoord();
    static TrajectoryIO* Alloc() { return (TrajectoryIO*)new Traj_AmberCoord(); }
  private:
    int natom3_;          ///< Number of coords (# atoms X 3)
    // TODO: Replace hasREMD with REMD_HEADER_SIZE/HasT()
    size_t hasREMD_;      ///< REMD header size if present
    int numBoxCoords_;    ///< Number of box coords, 3 (ortho or truncoct) or 6 (triclinic)
    const char* outfmt_;  ///< Format string for writing coordinates
    bool highPrecision_;  ///< If true output format will be 8.6 instead of 8.3
    BufferedFrame file_;  ///< Buffer reads/writes
    double boxAngle_[3];  ///< Hold default box angles in case traj has only box lengths

    static const size_t REMD_HEADER_SIZE;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void Info();
    int processWriteArgs(ArgList&);
    int readVelocity(int, double*);
    int processReadArgs(ArgList&)  { return 0; }
    int readIndices(int,int*)      { return 1; }
};
#endif
