#ifndef INC_TRAJECTORYFILEBASE_H
#define INC_TRAJECTORYFILEBASE_H
/// Class: TrajectoryFileBase
/// Base class that all trajectory types will inherit.
/// Currently TrajectoryFile and RemdTraj
#include "AmberParm.h" // BoxType.h
#include "ArgList.h"
#include "ProgressBar.h"
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
    int Frames;
    int numFramesRead;
    int total_read_frames;
    BoxType boxType;

public:
    TrajectoryFileBase();
    ~TrajectoryFileBase();
    void SetDebug(int);
    bool TrajFilenameIs(char *);
    // Inherited functions
    virtual int SetupRead(char *, ArgList *, AmberParm *) { return 1; }
    virtual int SetupWrite(char *, ArgList *, AmberParm *) { return 1; }
    virtual int BeginTraj(bool) { return 1; }
    virtual int EndTraj() { return 1; }
    virtual int GetNextFrame(double*,double*,double*) { return 1; }
    virtual int WriteFrame(int, AmberParm *, double*,double*,double) { return 1 };
    virtual void PrintInfo(int) { return; }
    // Functions that return private vars
    char *TrajName() { return trajName;}
    AmberParm *TrajParm() { return trajParm;}
    int Start() { return start; }
    int Total_Read_Frames() { return total_read_frames; }
    int NumFramesRead() { return numFramesRead; }
};
#endif
