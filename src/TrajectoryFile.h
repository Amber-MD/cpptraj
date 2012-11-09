#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
#include "TrajectoryIO.h"
/// Base class that all input and output trajectories will inherit.
class TrajectoryFile {
  public:
    /// Known trajectory formats.
    enum TrajFormatType {
      UNKNOWN_TRAJ=0, AMBERNETCDF, AMBERRESTARTNC, PDBFILE, MOL2FILE, CHARMMDCD,
      BINPOS, AMBERRESTART, AMBERTRAJ, CONFLIB
    };

    TrajectoryFile();
    virtual ~TrajectoryFile() {}

    /// Return string corresponding to given format
    static const char* FormatString( TrajFormatType );
    /// Get format type from keyword
    static TrajFormatType GetFormatFromArg(ArgList&);
    static TrajFormatType GetFormatFromString(std::string const&);
    /// Get standard file extension for trajectory format
    static std::string GetExtensionForType(TrajFormatType);
    /// Get type from extension
    static TrajFormatType GetTypeFromExtension(std::string const&);

    void SetDebug(int);
    void SetTrajFileName( std::string const& );
    int SetTrajParm( Topology* );

    Topology* TrajParm()              { return trajParm_;                  }
    const FileName& TrajName()        { return trajName_;                  }
    const char* BaseTrajStr()         { return trajName_.Base().c_str();   }
    const char* FullTrajStr()         { return trajName_.Full().c_str();   }
  protected:
    int debug_;            ///< Trajectory debug level.
    static TrajectoryIO* AllocTrajIO(TrajFormatType);
    static TrajectoryIO* DetectFormat(std::string const&);
    static TrajFormatType TrajFormat(std::string const&);
  private:
    struct TrajToken {
      TrajFormatType Type;
      const char* Key;
      const char* Description;
      const char* Extension;
      TrajectoryIO::AllocatorType Alloc;
    };
    static const TrajToken TrajArray[];
    typedef const TrajToken* TokenPtr;
    Topology *trajParm_;   ///< Associated parm
    FileName trajName_;    ///< The full path to trajectory file.
};
#endif
