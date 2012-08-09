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
      UNKNOWN_TRAJ=0, AMBERNETCDF, AMBERRESTARTNC, PDBFILE, MOL2FILE, CHARMMDCD,
      AMBERRESTART, AMBERTRAJ, CONFLIB, NTRAJ 
    };

    TrajectoryFile();
    ~TrajectoryFile();
    // Trajectory Setup functions
    int SetupRead(const char *, ArgList *, Topology *);
    int SetupWriteWithArgs(const char *, const char *, Topology *, TrajFormatType);
    // TODO: Accept string instead of char*
    int SetupWrite(const char *, ArgList *, Topology *, TrajFormatType);
    // Trajectory Read/Write functions
    int BeginTraj(bool);
    int EndTraj();
    int GetNextFrame(Frame&);
    int WriteFrame(int, Topology *, Frame&);

    void SetDebug(int);
    void SingleFrame();
    void PrintInfoLine();
    void PrintInfo(int);
    /// Get format type from keyword
    static TrajFormatType GetFormatFromArg(ArgList&);
    /// Return string corresponding to given format
    static const char* FormatString( TrajFormatType tIn );
    /// Get standard file extension for trajectory format
    static std::string GetExtensionForType(TrajFormatType);
    /// Get type from extension
    TrajFormatType GetTypeFromExtension(std::string const&);
    // Functions that return private vars
    const char *TrajName();
    Topology *TrajParm();
    int Start();
    int Total_Read_Frames();
    int Total_Frames();
    int NumFramesProcessed();
    bool HasVelocity();
    std::string FileName();
    const char* c_str();
  private:
    static const char FORMAT_STRINGS[][17];
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
    /// If defined, list of frame numbers to write.
    Range FrameRange_;
    /// If frame range define, this is next frame in range.
    Range::const_iterator rangeframe_;
    /// If true a frame range is defined.
    bool hasRange_;
    /// If true do not put box information in output traj
    bool nobox_;
    /// If true trajectory has been opened.
    bool trajIsOpen_;
    /// Set up TrajectoryIO object for reading multiple trajectories at once
    TrajectoryIO *setupRemdTrajIO(double, char*, TrajFormatType, ArgList&);
    /// Set up Trajectory IO object
    TrajectoryIO *SetupTrajectoryIO(TrajFormatType);
    /// Set up TrajectoryIO object for the given filename
    TrajectoryIO *setupTrajIO(const char *, TrajAccessType, TrajFormatType);
    /// Set start/stop/offset args from user input
    int SetArgs(ArgList *);
    /// Set actual start and stop
    int setupFrameInfo();
};
#endif
