#ifndef INC_PARM_MOL2_H
#define INC_PARM_MOL2_H
#include "ParmIO.h"
class Parm_Mol2 : public ParmIO {
  public :
    Parm_Mol2() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Mol2(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&) { return 1; }
    int processWriteArgs(ArgList&) { return 0; }
};
#endif
