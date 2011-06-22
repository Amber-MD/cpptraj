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
    enum BoxType { NONE, ORTHO, NONORTHO };

    int debug;
    TrajectoryIO *trajio;
    ProgressBar *progress;
    char *trajName;
    AmberParm *trajParm;
    AccessType fileAccess;
    int start;
    int stop;
    int offset;
    int Frames;
    BoxType boxType;
    // Specific to output traj
    Range *FrameRange;
    bool nobox;

    int setupTraj(char *, AccessType, FileFormat, FileType);
    void SetArgs(int, int, int);
  public:
    TrajectoryFile();
    ~TrajectoryFile();

    int SetupRead(char *, ArgList *, AmberParm *);
    int SetupWrite(char *, ArgList *);
    int BeginTraj();
    int EndTraj();
    int GetNextFrame();
    int WriteFrame(int, AmberParm *);
};
#endif
