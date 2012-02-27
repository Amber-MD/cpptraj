#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
class CharmmPsfParmFile : public ParmIO {
  public :
    int ReadParm(AmberParm&, CpptrajFile&);
};
#endif
