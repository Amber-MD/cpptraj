#ifndef INC_ACTION_DIHEDRAL_H
#define INC_ACTION_DIHEDRAL_H

#include "Action.h"

class Dihedral: public Action {
    DataSet *dih;
    AtomMask M1, M2, M3, M4;
  public:
    Dihedral();
    ~Dihedral();

    int init();
    int setup();
    int action();
};
#endif  
