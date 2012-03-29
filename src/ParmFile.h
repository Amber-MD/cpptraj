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

    void SetDebug(int);
    int Read(Topology&, char*,bool,bool);
    int Write(Topology&, char*,ParmFormatType);
  private :
    int debug_;
    ParmIO *SetupParmIO(ParmFormatType);
};
#endif
