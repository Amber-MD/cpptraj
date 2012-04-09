#ifndef INC_TRAJ_AMBERRESTART_H
#define INC_TRAJ_AMBERRESTART_H
#include "TrajectoryIO.h"
#include "FrameBuffer.h"
// Class: AmberRestart.h
/// Reads and writes formatted (ASCII text) amber
class AmberRestart : public TrajectoryIO, FrameBuffer {
  public:

    AmberRestart();
    //~AmberRestart();
    // AmberRestart-specific functions
    void SetNoVelocity();
  private:
    static const size_t BUF_SIZE;

    int restartAtoms_;     ///< Number of atoms in restart file
    int natom3_;           ///< Number of coords
    //size_t frameSize_;        ///< Size of 1 coord frame in bytes, inc box & velo if present
    //char *frameBuffer_;    ///< Used to read in restart coord
    int numBoxCoords_;     ///< Number of box coords (3 or 6)
    double restartTime_;   ///< Time in restart file, read in
    double restartTemp_;   ///< (Optional) replica temperature, read in.
    double time0_;         ///< For writes, restart time offset
    double dt_;            ///< For writes, restart timestep (scaling)
    bool singleWrite_;     ///< If false, frame # will be appended to output filename

    // Inherited functions
    bool ID_TrajFormat();
    int setupTrajin(Topology*);
    int setupTrajout(Topology*);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    int processWriteArgs(ArgList*);
    void info();

    int getBoxAngles(char *);
};
#endif
