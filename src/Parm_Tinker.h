#ifndef INC_PARM_TINKER_H
#define INC_PARM_TINKER_H
#include "ParmIO.h"
class Parm_Tinker : public ParmIO {
  public :
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Tinker(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&) { return 1; }
    int processWriteArgs(ArgList&) { return 0; }
};
#endif
