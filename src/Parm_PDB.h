#ifndef INC_PARM_PDB_H
#define INC_PARM_PDB_H
#include "ParmIO.h"
class PdbParmFile : public ParmIO {
  public :
    int ReadParm(AmberParm&, CpptrajFile&);
};
#endif
