#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    static ParmIO* Alloc() { return (ParmIO*)new Parm_CharmmPsf(); }
    bool ID_ParmFormat(CpptrajFile&);
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1; }
    void SetDebug(int) {}
};
#endif
