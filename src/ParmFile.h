#ifndef INC_PARMFILE_H
#define INC_PARMFILE_H
#include "ParmIO.h"
#include "Topology.h"
class ParmFile {
  public :
    enum ParmFormatType {
      UNKNOWN_PARM=0, PDBFILE, AMBERPARM, MOL2FILE, CHARMMPSF, NPARM
    };

    ParmFile();

    inline std::string &BaseName() {
      return basename_;
    }

    void SetDebug(int);
    int Read(Topology&, char*,bool,bool);
    int Write(Topology&, char*,ParmFormatType);
    int Write(Topology&, std::string, ParmFormatType);
  private :
    int debug_;
    std::string basename_;
    ParmIO *SetupParmIO(ParmFormatType);
};
#endif
