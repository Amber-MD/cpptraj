#ifndef INC_ACTION_DIHEDRAL_H
#define INC_ACTION_DIHEDRAL_H
#include "Action.h"
/// Calculate dihedral in a Frame
class Dihedral: public Action {
  public:
    Dihedral();
  private:
    DataSet *dih;
    AtomMask M1, M2, M3, M4;

    int init();
    int setup();
    int action();
};
#endif  
