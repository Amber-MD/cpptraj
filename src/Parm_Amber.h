#ifndef INC_PARM_AMBER_H
#define INC_PARM_AMBER_H
#include "ParmIO.h"
class AmberParmFile : public ParmIO {
  public :
    int ReadParm(AmberParm&, CpptrajFile&);
    int WriteParm(AmberParm&, char*);
  private :
    void SetParmFromValues(AmberParm &, int *, bool);
    int ReadParmAmber(AmberParm&, CpptrajFile&);
    int ReadParmOldAmber(AmberParm&, CpptrajFile&);
};
#endif
