#ifndef INC_PARM_PDB_H
#define INC_PARM_PDB_H
#include "ParmIO.h"
class Parm_PDB : public ParmIO {
  public :
    Parm_PDB() { }
    static ParmIO* Alloc() { return (ParmIO*)new Parm_PDB(); }
    bool ID_ParmFormat(CpptrajFile&);
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1; }
    void SetDebug(int) {}
};
#endif

