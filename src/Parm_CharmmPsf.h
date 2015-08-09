#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    Parm_CharmmPsf() : debug_(0) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CharmmPsf(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&);
    void SetDebug(int d) { debug_ = d; }
    int processWriteArgs(ArgList&) { return 0; }
    bool NeedsBondSearch() const { return false; }
  private:
    static inline int FindTag(char*, const char*, int, CpptrajFile&); 
    int debug_;
};
#endif
