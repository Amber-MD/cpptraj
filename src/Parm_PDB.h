#ifndef INC_PARM_PDB_H
#define INC_PARM_PDB_H
#include "ParmIO.h"
#include "PDBtype.h"
class Parm_PDB : public ParmIO, PDBtype {
  public :
    Parm_PDB() { }
    int ReadParm(Topology&);
    bool ID_ParmFormat();
};
#endif

