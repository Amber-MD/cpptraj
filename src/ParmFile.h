#ifndef INC_PARMFILE_H
#define INC_PARMFILE_H
#include "ParmIO.h"
#include "FileTypes.h"
class ParmFile {
    /// Allocator and description for file types. 
    static const FileTypes::AllocToken PF_AllocArray[];
    /// For associating keywords/extensions with file types. 
    static const FileTypes::KeyToken PF_KeyArray[];
    /// Topologies for which writes are supported.
    static const FileTypes::KeyToken PF_WriteKeyArray[];
  public :
    enum ParmFormatType { AMBERPARM=0, PDBFILE, MOL2FILE, CHARMMPSF, CIFFILE,
                          GMXTOP, SDFFILE, TINKER, UNKNOWN_PARM };
    static void ReadOptions(std::string const& fkey) { FileTypes::Options(PF_KeyArray,PF_AllocArray,UNKNOWN_PARM,fkey,FileTypes::READOPT); }
    static void WriteOptions(std::string const& fkey){ FileTypes::Options(PF_WriteKeyArray,PF_AllocArray,UNKNOWN_PARM,fkey,FileTypes::WRITEOPT);}
    ParmFile() {}
    int ReadTopology(Topology&, FileName const&, ArgList const&,int);
    int ReadTopology(Topology& t, FileName const& n, int d) {
      return ReadTopology(t, n, ArgList(), d);
    }
    int WritePrefixTopology(Topology const&, std::string const&, ParmFormatType,int);
    int WriteTopology(Topology const&, FileName const&, ArgList const&,ParmFormatType,int);
    int WriteTopology(Topology const& t, FileName const& n, ParmFormatType f,int d) {
      return WriteTopology(t, n, ArgList(), f, d);
    }
    FileName const& ParmFilename() { return parmName_; }
    /// \return ParmFormatType of given file or UNKNOWN_PARM.
    static ParmFormatType DetectFormat(FileName const&);
  private :
    /// \return Allocated ParmIO if given file matches known type, 0 otherwise.
    static ParmIO* DetectFormat(FileName const&, ParmFormatType&);

    FileName parmName_; ///< Topology input/output file name. 
};
#endif
