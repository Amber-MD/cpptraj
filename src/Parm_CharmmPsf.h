#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
#include "DataSet_Parameters.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    Parm_CharmmPsf() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CharmmPsf(); }
    bool ID_ParmFormat(CpptrajFile&);
    static void ReadHelp();
    int processReadArgs(ArgList&);
    int ReadParm(FileName const&, Topology&);
    int WriteParm(FileName const&, Topology const&);
    int processWriteArgs(ArgList&) { return 0; }
  private:
    static inline int FindTag(char*, const char*, int, CpptrajFile&);
    int ReadDihedrals(CpptrajFile&, int, const char*, Topology&) const;

    DataSet_Parameters params_;
};
#endif
