#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
/// Class: TrajectoryFile
/// Main trajectory class. Will have a TrajectoryIO class that will interface
/// with actual trajectory format.
#include "AmberParm.h" // BoxType.h
#include "ArgList.h"
#include "ProgressBar.h"
#include "TrajectoryIO.h"
#include "Range.h"
class TrajectoryFile {
    int debug;
    ProgressBar *progress;
    TrajectoryIO *trajio;
    char *trajName;
    AmberParm *trajParm;
    AccessType fileAccess;
    int start;
    int stop;
    int offset;
    int total_frames;      // Total # of frames in the traj
    int numFramesRead;     // # frames actually read
    int total_read_frames; // # frames that will be read based on start/stop/offset
    BoxType boxType;
    int currentFrame;
    // Specific to input traj
    int targetSet;
    int frameskip;
    // Specific to output traj
    Range *FrameRange;
    bool nobox;
    bool setupForWrite;

    TrajectoryIO *setupRemdTraj(char *, ArgList*);
    TrajectoryIO *setupTrajIO(char *, AccessType, FileFormat, FileType);
    void SetArgs(int, int, int);
    void SetTrajName(char *);
    int SetBoxType(TrajectoryIO *);

  public:
    TrajectoryFile();
    ~TrajectoryFile();
    // Trajectory IO functions
    int SetupRead(char *, ArgList *, AmberParm *);
    int SetupWrite(char *, ArgList *, AmberParm *);
    int BeginTraj(bool);
    int EndTraj();
    int GetNextFrame(double*,double*,double*);
    int WriteFrame(int, AmberParm *, double*,double*,double);
    // Public functions
    void SetDebug(int);
    void SingleFrame();
    bool TrajFilenameIs(char *);
    void PrintInfo(int);
    int SetupFrameInfo();
    // Functions that return private vars
    int CurrentFrame() { return currentFrame; }
    char *TrajName() { return trajName;}
    AmberParm *TrajParm() { return trajParm;}
    int Start() { return start; }
    int Total_Read_Frames() { return total_read_frames; }
    int Total_Frames() { return total_frames; }
    int NumFramesRead() { return numFramesRead; }
};
#endif
