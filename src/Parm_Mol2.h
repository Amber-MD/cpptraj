#ifndef INC_PARM_MOL2_H
#define INC_PARM_MOL2_H
#include "ParmIO.h"
class Parm_Mol2 : public ParmIO {
  public :
    static ParmIO* Alloc() { return (ParmIO*)new Parm_Mol2(); }
    bool ID_ParmFormat(CpptrajFile&);
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1; }
    void SetDebug(int) {}
};
#endif
