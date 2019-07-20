#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    Parm_CharmmPsf() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CharmmPsf(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&);
    int processWriteArgs(ArgList&) { return 0; }
  private:
    static const unsigned int ChmStrMax_;
    static inline int FindTag(char*, const char*, int, CpptrajFile&);
    static inline int ParseResID(char&, const char*);
};
#endif
