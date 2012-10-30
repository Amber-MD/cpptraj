#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
#include "TrajectoryIO.h"
/// Base class that all input and output trajectories will inherit.
class TrajectoryFile {
  public:
    /// Known trajectory formats.
    // NOTE: FORMAT_STRINGS must also be updated
    enum TrajFormatType {
      UNKNOWN_TRAJ=0, AMBERNETCDF, AMBERRESTARTNC, PDBFILE, MOL2FILE, CHARMMDCD,
      BINPOS, AMBERRESTART, AMBERTRAJ, CONFLIB, NTRAJ
    };

    TrajectoryFile();
    virtual ~TrajectoryFile() {}

    /// Return string corresponding to given format
    static const char* FormatString( TrajFormatType );
    /// Get format type from keyword
    static TrajFormatType GetFormatFromArg(ArgList&);
    /// Get standard file extension for trajectory format
    static std::string GetExtensionForType(TrajFormatType);
    /// Get type from extension
    static TrajFormatType GetTypeFromExtension(std::string const&);

    void SetDebug(int);
    void SetFileNames( std::string const&, std::string const& );
    int SetTrajParm( Topology* );

    Topology* TrajParm()              { return trajParm_;           }
    const char* FullTrajStr()         { return trajName_.c_str();   }
    std::string const& FullTrajName() { return trajName_;           }
    const char* BaseTrajStr()         { return baseName_.c_str();   }
    std::string const& BaseTrajName() { return baseName_;           }
  protected:
    int debug_;            ///< Trajectory debug level.

    TrajectoryIO* AllocTrajIO(TrajFormatType);
    TrajectoryIO* DetectFormat(CpptrajFile&);
  private:
    /// Strings describing each TrajFormatType
    static const char FORMAT_STRINGS[][17];

    Topology *trajParm_;   ///< Associated parm
    std::string trajName_; ///< The full path to trajectory file.
    std::string baseName_; ///< The base trajectory file name.
};
#endif
