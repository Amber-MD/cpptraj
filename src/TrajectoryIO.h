#ifndef INC_TRAJECTORYIO_H
#define INC_TRAJECTORYIO_H
/// Class: TrajectoryIO
/// TrajectoryIO will be the base class for performing trajectory reading
/// and writing that all formats will inherit.
#include "PtrajFile.h" 
class TrajectoryIO {
    PtrajFile *tfile;
    char *title;
  public:
    TrajectoryIO();
    virtual ~TrajectoryIO(); // virtual since this class is inherited.

    // Inherited classes
    virtual int setupRead(int) { return 1; }
    virtual int setupWrite(int) { return 1; }
    virtual int openTraj() { return 1; }
    virtual int readFrame(int,double*,double*,double*) { return 1; }
    virtual int writeFrame(int,double*,double*,double*) { return 1; }
    virtual int closeTraj() { return 1; }

    void SetFile(PtrajFile *);
    void SetTitle(char *);
}; 
#endif
