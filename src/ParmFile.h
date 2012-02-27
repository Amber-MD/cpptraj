#ifndef INC_PARMFILE_H
#define INC_PARMFILE_H
#include "AmberParm.h"
class ParmFile {
  public :
    ParmFile();
    ParmFile(int);

    int Read(AmberParm&, char*,bool,bool);
  private :
    int debug;
};
#endif
