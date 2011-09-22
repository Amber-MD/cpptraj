#ifndef INC_TRAJ_CHARMMDCD_H
#define INC_TRAJ_CHARMMDCD_H
#include "TrajectoryIO.h"
/// Class: CharmmDcd
/// TrajectoryIO class for reading coordinates from charmm dcd files.
class CharmmDcd : public TrajectoryIO {
    int dcdatom;

    union byte {
      unsigned char c[4];
      int i;
      float f;
    };

    // Inherited functions
    int setupRead(AmberParm *);
    int setupWrite(AmberParm *);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList *);

    int binaryByteOrder(int, int *);
    int readBinaryInteger(int);
  public :
    CharmmDcd();
    ~CharmmDcd();
    // charmm dcd-specific functions
};
#endif
