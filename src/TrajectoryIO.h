#ifndef INC_TRAJECTORYIO_H
#define INC_TRAJECTORYIO_H
/// Class: TrajectoryIO
/// TrajectoryIO will be the base class for performing trajectory reading
/// and writing that all formats will inherit.
#include "PtrajFile.h" 
class TrajectoryIO {
  protected:
    PtrajFile *tfile;
    char *title;
    bool hasTemperature;
    int debug;
  public:
    bool seekable;
    bool hasBox;
    double boxAngle[3]; // Hold alpha, beta and gamma angles

    TrajectoryIO();
    virtual ~TrajectoryIO(); // virtual since this class is inherited.

    // Inherited functions 
    virtual int setupRead(int) { return 1; }
    virtual int setupWrite(int) { return 1; }
    virtual int openTraj() { return 1; }
    virtual int readFrame(int,double*,double*,double*) { return 1; }
    virtual int writeFrame(int,double*,double*,double) { return 1; }
    virtual void closeTraj() { return; }
    virtual void info() { return; }

    void SetFile(PtrajFile *);
    void SetTitle(char *);
    bool FilenameIs(char *); 
}; 
#endif
