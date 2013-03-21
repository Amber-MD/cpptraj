#ifndef INC_PARMFILE_H
#define INC_PARMFILE_H
#include "ParmIO.h"
class ParmFile {
  public :
    enum ParmFormatType { UNKNOWN_PARM=0, PDBFILE, AMBERPARM, MOL2FILE, CHARMMPSF };

    ParmFile() {}
    int Read(Topology&, std::string const&, bool,int);
    int Write(Topology const&, std::string const&, ParmFormatType,int);
    FileName const& ParmFilename() { return parmName_; }
  private :
    struct ParmToken {
      ParmFormatType Type;
      const char* Key;
      ParmIO::AllocatorType Alloc;
    };
    static const ParmToken ParmArray[];
    typedef const ParmToken* TokenPtr;
    FileName parmName_; ///< Used on Read for TopologyList
};
#endif
