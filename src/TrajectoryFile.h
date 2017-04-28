#ifndef INC_TRAJECTORYFILE_H
#define INC_TRAJECTORYFILE_H
#include "TrajectoryIO.h"
#include "FileTypes.h"
/// Base class that all input and output trajectories will inherit.
/** There are 3 steps to adding new trajectory types:
  *   - 1) Create the TrajectoryIO-derived class for the format and include
  *        it in TrajectoryFile.cpp.
  *   - 2) Add a unique entry to enumerated type TrajFormatType.
  *   - 3) Add entry/entries describing how the format is to be called
  *        to the TF_AllocArray[] and TF_KeyArray[] arrays.
  */
// TODO: Namespace? Combine detection of all format types to one file?
class TrajectoryFile {
    /// Allocator and description for file types. 
    static const FileTypes::AllocToken TF_AllocArray[];
    /// For associating keywords/extensions with file types. 
    static const FileTypes::KeyToken TF_KeyArray[];
  public:
    /// Known trajectory formats.
    enum TrajFormatType {
      AMBERNETCDF = 0, AMBERRESTARTNC, AMBERNCENSEMBLE, PDBFILE, MOL2FILE, CIF, CHARMMDCD, 
      GMXTRX, GMXXTC, BINPOS, AMBERRESTART, GRO, TINKER, CHARMMCOR, AMBERTRAJ, SQM, SDF, CONFLIB,
      UNKNOWN_TRAJ
    };

    TrajectoryFile() {}
    virtual ~TrajectoryFile() {}
    /// List read options for each format.
    static void ReadOptions() { FileTypes::ReadOptions(TF_KeyArray,TF_AllocArray, UNKNOWN_TRAJ); }
    /// List write options for each format.
    static void WriteOptions(){ FileTypes::WriteOptions(TF_KeyArray,TF_AllocArray,UNKNOWN_TRAJ); }
    /// \return format type from keyword in ArgList. 
    static TrajFormatType GetFormatFromArg(ArgList& a) {
      return (TrajFormatType)FileTypes::GetFormatFromArg(TF_KeyArray, a, UNKNOWN_TRAJ);
    }
    /// \return format type from keyword. Default to Amber trajectory.
    static TrajFormatType GetFormatFromString(std::string const& s) {
      return (TrajFormatType)FileTypes::GetFormatFromString(TF_KeyArray, s, AMBERTRAJ);
    }
    /// \return format type from keyword, or given default if keyword is empty.
    static TrajFormatType GetFormatFromString(std::string const& s, TrajFormatType def) {
      return (TrajFormatType)FileTypes::GetFormatFromString(TF_KeyArray, s, def);
    }
    /// \return standard file extension for trajectory format.
    static std::string GetExtensionForType(TrajFormatType t) {
      return FileTypes::GetExtensionForType(TF_KeyArray, t);
    }
    /// \return type from extension.
    static TrajFormatType GetTypeFromExtension(std::string const& e) {
      return (TrajFormatType)FileTypes::GetTypeFromExtension(TF_KeyArray, e, UNKNOWN_TRAJ);
    }
    /// \return type from extension. Default to given format if unknown.
    static TrajFormatType GetTypeFromExtension(FileName const& f, TrajFormatType def) {
      return (TrajFormatType)FileTypes::GetTypeFromExtension(TF_KeyArray, f.Ext(), def);
    }
    /// \return string corresponding to given format.
    static const char* FormatString( TrajFormatType tt ) { return 
      FileTypes::FormatDescription(TF_AllocArray, tt);
    }
    /// Allocate TrajectoryIO appropriate for given file.
    static TrajectoryIO* DetectFormat(FileName const&, TrajFormatType&);
    /// Allocate TrajectoryIO for given format
    static TrajectoryIO* AllocTrajIO(TrajFormatType t) {
      return (TrajectoryIO*)FileTypes::AllocIO(TF_AllocArray, t, true);
    }
};
#endif
