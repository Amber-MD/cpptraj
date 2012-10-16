#ifndef INC_ACTION_DIHEDRAL_H
#define INC_ACTION_DIHEDRAL_H
#include "Action.h"
/// Calculate dihedral in a Frame
class Action_Dihedral: public Action {
  public:
    Action_Dihedral();
  private:
    DataSet* dih_;
    bool useMass_;
    AtomMask M1_;
    AtomMask M2_;
    AtomMask M3_;
    AtomMask M4_;

    int init();
    int setup();
    int action();
};
#endif  
