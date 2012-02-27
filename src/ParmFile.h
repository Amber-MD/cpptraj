#ifndef INC_PARMFILE_H
#define INC_PARMFILE_H
#include "FileRoutines.h" // FileFormat
#include "AmberParm.h"
class ParmFile {
  public :
    ParmFile();

    void SetDebug(int);
    int Read(AmberParm&, char*,bool,bool);
    int Write(AmberParm&, char*,FileFormat);
  private :
    int debug;
};
#endif
