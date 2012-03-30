#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    bool ID_ParmFormat();
    int ReadParm(Topology&);
  private:
    static const size_t BUF_SIZE_ = 256;
    char buffer_[BUF_SIZE_];
};
#endif
