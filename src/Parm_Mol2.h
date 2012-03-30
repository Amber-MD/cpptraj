#ifndef INC_PARM_MOL2_H
#define INC_PARM_MOL2_H
#include "ParmIO.h"
#include "Mol2File.h"
class Parm_Mol2 : public ParmIO, Mol2File {
  public :
    bool ID_ParmFormat();
    int ReadParm(Topology&);
};
#endif
