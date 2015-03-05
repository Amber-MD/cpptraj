#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    Parm_CharmmPsf() : debug_(0) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CharmmPsf(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&);
    void SetDebug(int d) { debug_ = d; }
    int processWriteArgs(ArgList&) { return 0; }
  private:
    static inline int FindTag(char*, const char*, int, CpptrajFile&); 
    int debug_;
};
#endif
