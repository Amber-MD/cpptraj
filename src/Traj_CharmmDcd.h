#ifndef INC_TRAJ_CHARMMDCD_H
#define INC_TRAJ_CHARMMDCD_H
#include "TrajectoryIO.h"
/// Class: CharmmDcd
/// TrajectoryIO class for reading coordinates from charmm dcd files.
class CharmmDcd : public TrajectoryIO {
    int dcdatom;
    bool dcdExtraBlock;
    bool dcd4D;
    int istart;
    int nsavc;
    int namnf;
    int nfreat;
    int *freeat;
    float timestep;

    union doublebyte {
      unsigned char c[8];
      int i[2];
      double d;
    };
    union headerbyte {
      unsigned char c[80];
      int i[20];
      float f[20];
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

  public :
    CharmmDcd();
    ~CharmmDcd();
    // charmm dcd-specific functions
};
#endif
