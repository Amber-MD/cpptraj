#ifndef INC_ACTION_DIHEDRAL_H
#define INC_ACTION_DIHEDRAL_H
#include "Action.h"
/// Calculate dihedral in a Frame
class Dihedral: public Action {
    DataSet *dih;
    AtomMask M1, M2, M3, M4;
  public:
    Dihedral();

    int init();
    int setup();
    int action();
};
#endif  
