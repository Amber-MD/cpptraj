#ifndef INC_PARM_MOL2_H
#define INC_PARM_MOL2_H
#include "ParmIO.h"
class Parm_Mol2 : public ParmIO {
  public :
    Parm_Mol2() : needsBondSearch_(false) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Mol2(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1; }
    void SetDebug(int) {}
    int processWriteArgs(ArgList&) { return 0; }
    bool NeedsBondSearch() const { return needsBondSearch_; }
  private:
    bool needsBondSearch_;
};
#endif
