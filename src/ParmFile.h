#ifndef INC_PARMFILE_H
#define INC_PARMFILE_H
#include "ParmIO.h"
#include "AmberParm.h"
class ParmFile {
  public :
    enum ParmFormatType {
      UNKNOWN_PARM=0, PDBFILE, AMBERPARM, MOL2FILE, CHARMMPSF, OLDAMBERPARM
    };

    ParmFile();

    void SetDebug(int);
    int Read(AmberParm&, char*,bool,bool);
    int Write(AmberParm&, char*,ParmFormatType);
  private :
    int debug_;

    ParmFormatType ID_ParmFormat(ParmIO &);
};
#endif
