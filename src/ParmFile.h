#ifndef INC_PARMFILE_H
#define INC_PARMFILE_H
#include "FileTypes.h"
#include "FileName.h"
// Forward declares
class ParmIO;
class Topology;
/// Class that reads/writes Topology class.
class ParmFile {
    /// Allocator and description for file types. 
    static const FileTypes::AllocToken PF_AllocArray[];
    /// For associating keywords/extensions with file types. 
    static const FileTypes::KeyToken PF_KeyArray[];
    /// Topologies for which writes are supported.
    static const FileTypes::KeyToken PF_WriteKeyArray[];
  public :
    /// Recognized format types.
    enum ParmFormatType { AMBERPARM=0, PDBFILE, MOL2FILE, CHARMMPSF, CIFFILE,
                          GMXTOP, SDFFILE, TINKER, UNKNOWN_PARM };
    /// List of read options corresponding to format key.
    static void ReadOptions(std::string const& fkey) { FileTypes::Options(PF_KeyArray,PF_AllocArray,UNKNOWN_PARM,fkey,FileTypes::READOPT); }
    /// List of write options corresponding to format key.
    static void WriteOptions(std::string const& fkey){ FileTypes::Options(PF_WriteKeyArray,PF_AllocArray,UNKNOWN_PARM,fkey,FileTypes::WRITEOPT);}
    /// CONSTRUCTOR
    ParmFile();
    /// ReadTopology() keywords
    static const char* ReadTopologyKeywords();
    /// ReadTopology() help
    static const char* ReadTopologyHelp();
    /// Read topology file with optional arguments and debug level
    int ReadTopology(Topology&, FileName const&, ArgList const&,int);
    /// Read topology file
    int ReadTopology(Topology&, FileName const&, int);
    /// Write Topology to specified file as <prefix>.<originalFileName> with optional args
    int WritePrefixTopology(Topology const&, std::string const&, ArgList const&,ParmFormatType,int);
    /// Write Topology to specified file as <prefix>.<originalFileName>
    int WritePrefixTopology(Topology const&, std::string const&, ParmFormatType,int);
    /// Write Topology to specified file with optional arguments and debug level
    int WriteTopology(Topology const&, FileName const&, ArgList const&,ParmFormatType,int);
    /// Write Topology to specified file
    int WriteTopology(Topology const&, FileName const&, ParmFormatType, int);
    /// \return File name
    FileName const& ParmFilename() { return parmName_; }
    /// \return ParmFormatType of given file or UNKNOWN_PARM.
    static ParmFormatType DetectFormat(FileName const&);
  private :
    /// \return Allocated ParmIO if given file matches known type, 0 otherwise.
    static ParmIO* DetectFormat(FileName const&, ParmFormatType&);

    FileName parmName_; ///< Topology input/output file name.
};
#endif
