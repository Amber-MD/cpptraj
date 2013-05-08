#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
#include "TrajectoryIO.h"
/// Base class that all input and output trajectories will inherit.
/** There are 3 steps to adding new trajectory types:
  *   - 1) Create the TrajectoryIO-derived class for the format and include
  *        it in TrajectoryFile.cpp.
  *   - 2) Add a unique entry to enumerated type TrajFormatType.
  *   - 3) Add entry/entries describing how the format is to be called
  *        to the TrajArray[] array.
  */
class TrajectoryFile {
  public:
    /// Known trajectory formats.
    enum TrajFormatType {
      UNKNOWN_TRAJ=0, AMBERNETCDF, AMBERRESTARTNC, PDBFILE, MOL2FILE, CHARMMDCD,
      BINPOS, AMBERRESTART, AMBERTRAJ, CONFLIB, GMXTRX
    };

    TrajectoryFile();
    virtual ~TrajectoryFile() {}
    /// Return string corresponding to given format.
    static const char* FormatString( TrajFormatType );
    /// Get format type from keyword in ArgList. 
    static TrajFormatType GetFormatFromArg(ArgList&);
    /// Get format type from keyword.
    static TrajFormatType GetFormatFromString(std::string const&);
    /// Get standard file extension for trajectory format
    static std::string GetExtensionForType(TrajFormatType);
    /// Get type from extension
    static TrajFormatType GetTypeFromExtension(std::string const&);

    void SetDebug(int);
    void SetTrajFileName( std::string const&, bool );
    int SetTrajParm( Topology* );

    Topology* TrajParm()              { return trajParm_;                }
    const FileName& TrajFilename()    { return trajName_;                }
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
