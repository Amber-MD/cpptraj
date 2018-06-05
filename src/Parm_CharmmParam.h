#ifndef INC_PARM_CHARMMPARAM_H
#define INC_PARM_CHARMMPARAM_H
#include "ParmIO.h"
class Parm_CharmmParam : public ParmIO {
  public :
    Parm_CharmmParam() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CharmmParam(); }
    bool ID_ParmFormat(CpptrajFile&) { return false; }
    //static void ReadHelp();
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&) { return 1; }
    int WriteParm(FileName const&, Topology const&);
    int processWriteArgs(ArgList&) { return 0; }
  private:
};
#endif
