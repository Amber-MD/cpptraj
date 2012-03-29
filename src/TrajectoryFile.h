#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
#include "ProgressBar.h"
#include "TrajectoryIO.h" // Topology BoxType CpptrajFile ArgList 
#include "Range.h"
#include "Frame.h"
// Class: TrajectoryFile
/// Allow reading and writing of trajectory files.
/** Wrapper around the TrajectoryIO base class that allows the rest of cpptraj
  * to interface with trajectories. Holds general trajectory information such
  * as total # of frames, read start, stop, and offset values, etc.
  */
class TrajectoryFile {
  public:
    enum TrajFormatType {
      UNKNOWN_TRAJ=0, PDBFILE, AMBERTRAJ, AMBERNETCDF, AMBERRESTART,
      CONFLIB, AMBERRESTARTNC, MOL2FILE, CHARMMDCD
    };

    TrajectoryFile();
    ~TrajectoryFile();
    // Trajectory IO functions
    int SetupRead(char *, ArgList *, Topology *);
    int SetupWriteWithArgs(char *, const char *, Topology *, TrajFormatType);
    int SetupWrite(char *, Topology *, char *);
    // Next two currently only used for Clustering
    int SetupWrite(char *, ArgList *, Topology *, TrajFormatType);
    int SetupNumberedWrite(char *, int, Topology *, char *);
    int BeginTraj(bool);
    int EndTraj();
    int GetNextFrame(Frame&);
    int WriteFrame(int, Topology *, Frame&);
    // Public functions
    void SetDebug(int);
    void SingleFrame();
    bool TrajFilenameIs(char *);
    void PrintInfoLine();
    void PrintInfo(int);
    // Functions that return private vars
    int CurrentFrame();
    char *TrajName();
    Topology *TrajParm();
    int Start();
    int Total_Read_Frames();
    int Total_Frames();
    int NumFramesProcessed();
    bool HasVelocity();

    const char *FileName();
  private:
    /// Denote whether reading, writing, or appending.
    enum TrajAccessType { READTRAJ, WRITETRAJ, APPENDTRAJ };
    /// Trajectory debug level
    int debug_;
    /// Keep track of trajectory progress
    ProgressBar *progress_;
    /// Class that performs the actual IO for trajectory format
    TrajectoryIO *trajio_;
    /// Trajectory name
    const char *trajName_;
    /// Associated parm
    Topology *trajParm_;
    /// Trajectory File access (R/W/A)
    TrajAccessType fileAccess_;
    /// Number of frames that have been read/written with this traj
    int numFramesProcessed_;
    // -----===== Specific to input traj =====-----
    /// Frame to begin processing
    int start_;
    /// Frame to end processing
    int stop_;
    /// Number of frames to skip between processed frames
    int offset_;
    /// The total number of frames in the traj
    int total_frames_;
    /// The number of frames that will actually be read based on start/stop/offset
    int total_read_frames_;
    /// The current frame number being read
    int currentFrame_;
    /// The next frame to process
    int targetSet_;
    /// The number of frames to skip over while reading
    int frameskip_;
    // -----===== Specific to output traj =====-----
    /// If defined, list of frame numbers to write
    Range *FrameRange_;
    /// If true do not put box information in output traj
    bool nobox_;
    /// If true trajectory has been opened.
    bool trajIsOpen_;
    /// Identify trajectory format
    TrajFormatType ID_TrajFormat(TrajectoryIO &);
    /// Set up TrajectoryIO object for reading multiple trajectories at once
    TrajectoryIO *setupRemdTrajIO(double, char*, TrajFormatType, ArgList&);
    /// Set up TrajectoryIO object for the given filename
    TrajectoryIO *setupTrajIO(char *, TrajAccessType, TrajFormatType);
    /// Set start/stop/offset args from user input
    int SetArgs(ArgList *);
    /// Set actual start and stop
    int setupFrameInfo();
    /// Get format type from keyword
    TrajFormatType GetFormatFromArg(ArgList*);
    /// Get standard file extension for trajectory format
    std::string GetExtensionForType(TrajFormatType);
};
#endif
