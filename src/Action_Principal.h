#ifndef INC_ACTION_PRINCIPAL_H
#define INC_ACTION_PRINCIPAL_H
#include "Action.h"
#include "Lapack_Diag.h"
class Action_Principal : public Action {
  public:
    Action_Principal();
  private:
    Lapack_Diag Principal_;
    bool doRotation_;
    AtomMask mask_;

    int init();
    int setup();
    int action();
};
#endif
