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
    /// Known trajectory formats.
    // NOTE: FORMAT_STRINGS must also be updated
    enum TrajFormatType {
      UNKNOWN_TRAJ=0, AMBERNETCDF, AMBERRESTARTNC, PDBFILE, MOL2FILE, CHARMMDCD,
      BINPOS, AMBERRESTART, AMBERTRAJ, CONFLIB, NTRAJ 
    };

    TrajectoryFile();
    ~TrajectoryFile();
    // Trajectory Setup functions
    int SetupTrajRead(std::string const&, ArgList *, Topology *);
    int SetupTrajWriteWithArgs(std::string const&, const char *, Topology *, TrajFormatType);
    int SetupTrajWrite(std::string const&, ArgList *, Topology *, TrajFormatType);
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
    /// Indicate whether underlying TrajectoryIO object has velocity info.
    bool HasVelocity();
    // Functions that return private vars
    Topology* TrajParm()          { return trajParm_;           }
    int Start()                   { return start_;              }
    int Total_Read_Frames()       { return total_read_frames_;  }
    int Total_Frames()            { return total_frames_;       }
    int NumFramesProcessed()      { return numFramesProcessed_; }
    const char* FullTrajStr()         { return trajName_.c_str();   }
    std::string const& FullTrajName() { return trajName_;           }
    const char* BaseTrajStr()         { return baseName_.c_str();   }
    std::string const& BaseTrajName() { return baseName_;           }
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
    TrajectoryIO *velio_;
    /// The full path to trajectory file.
    std::string trajName_;
    /// The base trajectory file name.
    std::string baseName_;
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
    TrajectoryIO* setupRemdTrajIO(double, const char*, TrajFormatType, ArgList&, std::vector<int> const&);
    /// Set up Trajectory IO object
    TrajectoryIO* SetupTrajectoryIO(TrajFormatType);
    /// Set up TrajectoryIO object for the given filename
    TrajectoryIO* setupTrajIO(const char *, TrajAccessType, TrajFormatType);
    /// Set start/stop/offset args from user input
    int SetArgs(ArgList *);
    /// Set actual start and stop
    int setupFrameInfo();
};
#endif
