#ifndef INC_TRAJ_CHARMMDCD_H
#define INC_TRAJ_CHARMMDCD_H
#include "TrajectoryIO.h"
// Class: CharmmDcd
/// TrajectoryIO class for reading coordinates from charmm dcd files.
class CharmmDcd : public TrajectoryIO {
  public :
    CharmmDcd();
    ~CharmmDcd();
    // charmm dcd-specific functions
  private:
    int dcdatom;
    int dcdframes;
    int dcdoutsize;
    int dcdheadersize;
    bool isBigEndian;
    bool is64bit;
    unsigned int readSize;
    //bool dcdExtraBlock;
    bool dcd4D;
    int istart;
    int nsavc;
    int namnf;
    int nfreat;
    int *freeat;
    float timestep;
    float *xcoord;
    float *ycoord;
    float *zcoord;

    union headerbyte {
      unsigned char c[80];
      int i[20];
      float f[20];
    };
    int ReadBlock(int);
    int WriteBlock(int);
    int readDcdHeader();
    int writeDcdHeader();

    // Inherited functions
    bool ID_TrajFormat();
    int setupTrajin(Topology *);
    int setupTrajout(Topology *,int);
    int openTraj();
    void closeTraj();
    int readFrame(int,double*,double*,double*,double*);
    int writeFrame(int,double*,double*,double*,double);
    void info();
    int processWriteArgs(ArgList *);

};
#endif
