#ifndef INC_PARM_MOL2_H
#define INC_PARM_MOL2_H
#include "ParmIO.h"
class Mol2ParmFile : public ParmIO {
  public :
    int ReadParm(AmberParm&, CpptrajFile&);
};
#endif
