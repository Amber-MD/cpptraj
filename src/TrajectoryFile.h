#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
#include "ProgressBar.h"
#include "TrajectoryIO.h" // AmberParm BoxType CpptrajFile ArgList 
#include "Range.h"
#include "Frame.h"
// Class: TrajectoryFile
/// Allow reading and writing of trajectory files.
/** Wrapper around the TrajectoryIO base class that allows the rest of cpptraj
  * to interface with trajectories. Holds general trajectory information such
  * the start, stop, and offset values.
  */
// NOTE: Remove boxType?
class TrajectoryFile {
    /// trajectory debug level
    int debug;
    /// Keep track of trajectory progress
    ProgressBar *progress;
    /// Class that performs the actual IO for trajectory format
    TrajectoryIO *trajio;
    /// Trajectory name
    char *trajName;
    /// Associated parm
    AmberParm *trajParm;
    /// Trajectory File access (R/W/A)
    AccessType fileAccess;
    /// Number of frames that have been read/written with this traj
    int numFramesProcessed;
    /// True when the trajectory is open and ready for read/write 
    bool trajectoryIsOpen;
    // -----===== Specific to input traj =====-----
    /// Frame to begin processing
    int start;
    /// Frame to end processing
    int stop;
    /// Number of frames to skip between processed frames
    int offset;
    /// The total number of frames in the traj
    int total_frames;
    /// The number of frames that will actually be read based on start/stop/offset
    int total_read_frames;
    /// The box type of this trajectory if box coords present
    BoxType boxType;
    /// The current frame number being read
    int currentFrame;
    /// The next frame to process
    int targetSet;
    /// The number of frames to skip over while reading
    int frameskip;
    // -----===== Specific to output traj =====-----
    /// If defined, list of frame numbers to write
    Range *FrameRange;
    /// If true do not put box information in output traj
    bool nobox;
    /// Set up TrajectoryIO object for reading multiple trajectories at once
    TrajectoryIO *setupRemdTrajIO(char *, double, char*, FileFormat, ArgList&);
    /// Set up TrajectoryIO object for the given filename
    TrajectoryIO *setupTrajIO(char *, AccessType, FileFormat, FileType);
    /// Set start/stop/offset args from user input
    int SetArgs(ArgList *);
    /// Set the trajectory name
    void SetTrajName(char *);
    /// Set the trajectory box type
    int SetBoxType(TrajectoryIO *);
    int setupFrameInfo();
    FileFormat getFmtFromArg(ArgList*,FileFormat);

  public:
    TrajectoryFile();
    ~TrajectoryFile();
    // Trajectory IO functions
    int SetupRead(char *, ArgList *, AmberParm *);
    int SetupWriteWithArgs(char *, const char *, AmberParm *, FileFormat);
    int SetupWrite(char *, ArgList *, AmberParm *, FileFormat);
    int BeginTraj(bool);
    int EndTraj();
    int GetNextFrame(Frame&);
    int WriteFrame(int, AmberParm *, Frame&);
    // Public functions
    void SetDebug(int);
    void SingleFrame();
    bool TrajFilenameIs(char *);
    void PrintInfoLine();
    void PrintInfo(int);
    // Functions that return private vars
    int CurrentFrame() { return currentFrame; }
    char *TrajName() { return trajName;}
    AmberParm *TrajParm() { return trajParm;}
    int Start() { return start; }
    int Total_Read_Frames() { return total_read_frames; }
    int Total_Frames() { return total_frames; }
    int NumFramesProcessed() { return numFramesProcessed; }
    bool HasVelocity();
};
#endif
