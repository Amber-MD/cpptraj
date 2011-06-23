#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
/// Class: TrajectoryFile
/// Main trajectory class. Will have a TrajectoryIO class that will interface
/// with actual trajectory format.
#include "TrajectoryIO.h"
#include "AmberParm.h" // BoxType.h
#include "ArgList.h"
#include "ProgressBar.h"
#include "Range.h"
class TrajectoryFile {
    int debug;
    TrajectoryIO *trajio;
    FileFormat trajFormat;
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
    int currentFrame;
    // Specific to input traj
    int targetSet;
    int frameskip;
    // Specific to output traj
    Range *FrameRange;
    bool nobox;
    bool setupForWrite;

    int setupTraj(char *, AccessType, FileFormat, FileType);
    void SetArgs(int, int, int);
  public:
    TrajectoryFile();
    ~TrajectoryFile();

    int SetupRead(char *, ArgList *, AmberParm *);
    int SetupFrameInfo();
    int SetupWrite(char *, ArgList *, AmberParm *);
    int BeginTraj(bool);
    int EndTraj();
    int GetNextFrame(double*,double*,double*);
    int WriteFrame(int, AmberParm *, double*,double*,double);

    bool TrajFilenameIs(char *);
};
#endif
