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
    /// Trajectories for which writes are supported
    static const FileTypes::KeyToken TF_WriteKeyArray[];
  public:
    /// Known trajectory formats.
    enum TrajFormatType {
      AMBERNETCDF = 0, AMBERRESTARTNC, AMBERNCENSEMBLE, PDBFILE, MOL2FILE, CIF, CHARMMDCD, 
      GMXTRX, GMXXTC, BINPOS, AMBERRESTART, GRO, TINKER, CHARMMCOR, CHARMMREST,AMBERTRAJ,
      SQM, SDF, XYZ, CONFLIB,
      UNKNOWN_TRAJ
    };

    TrajectoryFile() {}
    virtual ~TrajectoryFile() {}
    /// List read options for each format.
    static void ReadOptions() { FileTypes::ReadOptions(TF_KeyArray,TF_AllocArray, UNKNOWN_TRAJ); }
    /// List write options for each format.
    static void WriteOptions(){ FileTypes::WriteOptions(TF_WriteKeyArray,TF_AllocArray,UNKNOWN_TRAJ); }
    /// \return write format type corresponding to given string, or default if no match.
    static TrajFormatType WriteFormatFromString(std::string const& s, TrajFormatType def) {
      return (TrajFormatType)FileTypes::GetFormatFromString(TF_WriteKeyArray,s,def);
    }
    /// \return write format type corresponding to extension of give filename, or default.
    static TrajFormatType WriteFormatFromFname(FileName const& f, TrajFormatType def) {
      return (TrajFormatType)FileTypes::GetTypeFromExtension(TF_WriteKeyArray,f.Ext(),def);
    }
    /// \return default filename extension for given write format type.
    static std::string WriteFormatExtension(TrajFormatType t) {
      return FileTypes::GetExtensionForType(TF_WriteKeyArray, t);
    }
    /// \return write format type corresponding to keyword in ArgList, or default type.
    static TrajFormatType WriteFormatFromArg(ArgList& a, TrajFormatType def) {
      return (TrajFormatType)FileTypes::GetFormatFromArg(TF_WriteKeyArray, a, def);
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
    /// \return TrajFormatType of given file or UNKNOWN_TRAJ
    static TrajFormatType DetectFormat(FileName const&);
    /// \return Allocated TrajectoryIO for given file with optional format keyword.
    static TrajectoryIO* DetectFormat(FileName const&, std::string const&, TrajFormatType&);
};
#endif
