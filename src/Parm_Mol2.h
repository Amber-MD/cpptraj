#ifndef INC_PARM_MOL2_H
#define INC_PARM_MOL2_H
#include "ParmIO.h"
#include "Mol2File.h"
class Mol2ParmFile : public ParmIO, Mol2File {
  public :
    int ReadParm(Topology&);
};
#endif
