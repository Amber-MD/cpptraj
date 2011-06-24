#ifndef INC_TRAJECTORYFILEBASE_H
#define INC_TRAJECTORYFILEBASE_H
/// Class: TrajectoryFileBase
/// Base class that all trajectory types will inherit.
/// Currently TrajectoryFile and RemdTraj
#include "AmberParm.h" // BoxType.h
#include "ArgList.h"
#include "ProgressBar.h"
#include "TrajectoryIO.h"
class TrajectoryFileBase {
  protected:
    int debug;
    ProgressBar *progress;
    char *trajName;
    AmberParm *trajParm;
    AccessType fileAccess;
    int start;
    int stop; 
    int offset;
    int total_frames;      // Total # of frames in the traj
    int numFramesRead;     // # frames actually read
    int total_read_frames; // # frames that will be read based on start/stop/offset
public:
    TrajectoryFileBase();
    virtual ~TrajectoryFileBase(); // virtual since class will be inherited
    void SetDebug(int);
    void SetTrajName(char *);
    void SetArgs(int, int, int);
    void SingleFrame();
    // Inherited functions
    virtual bool TrajFilenameIs(char *) {return false;}
    virtual int SetupRead(char *, ArgList *, AmberParm *) { return 1; }
    virtual int SetupWrite(char *, ArgList *, AmberParm *) { return 1; }
    virtual int BeginTraj(bool) { return 1; }
    virtual int EndTraj() { return 1; }
    virtual int GetNextFrame(double*,double*,double*) { return 1; }
    virtual int WriteFrame(int, AmberParm *, double*,double*,double) { return 1; }
    virtual void PrintInfo(int) { return; }
    virtual int SetupFrameInfo(); 
    // Functions that return private vars
    char *TrajName() { return trajName;}
    AmberParm *TrajParm() { return trajParm;}
    int Start() { return start; }
    int Total_Read_Frames() { return total_read_frames; }
    int Total_Frames() { return total_frames; }
    int NumFramesRead() { return numFramesRead; }
};
#endif
