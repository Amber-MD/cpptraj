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
                          GMXTOP, SDFFILE, TINKER, CHARMMPARM, UNKNOWN_PARM };
    static void ReadOptions() { FileTypes::ReadOptions(PF_KeyArray,PF_AllocArray,UNKNOWN_PARM); }
    static void WriteOptions(){ FileTypes::WriteOptions(PF_WriteKeyArray,PF_AllocArray,UNKNOWN_PARM);}
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
  private :
    ParmIO* DetectFormat(FileName const&, ParmFormatType&); 
    FileName parmName_; ///< Topology input/output file name. 
};
#endif
