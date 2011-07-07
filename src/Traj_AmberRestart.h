#ifndef INC_TRAJ_AMBERRESTART_H
#define INC_TRAJ_AMBERRESTART_H
#include "TrajectoryIO.h"
/// Class: AmberRestart.h
/// TrajectoryIO class for reading and writing formatted (ASCII text) amber
/// restart files.
class AmberRestart : public TrajectoryIO {
    int restartAtoms;     // Number of atoms in restart file
    int natom3;           // Number of coords
    int frameSize;        // Size of 1 coord frame in bytes, inc box & velo if present
    char *frameBuffer;    // Used to read in restart coord
    int numBoxCoords;     // Number of box coords (3 or 6)
    double restartTime;   // Time in restart file, read in
    double restartTemp;   // (Optional) replica temperature, read in.

    // Inherited functions
    int setupRead(AmberParm*);
    int setupWrite(AmberParm*);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();

    void getBoxAngles(char *, int);

  public:

  AmberRestart();
  ~AmberRestart();
  // AmberRestart-specific functions
};
#endif
