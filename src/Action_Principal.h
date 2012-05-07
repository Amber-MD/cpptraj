#ifndef INC_ACTION_PRINCIPAL_H
#define INC_ACTION_PRINCIPAL_H
#include "Action.h"
//#incl ude "Lapack_Diag.h"
//#incl ude "Lapack_General.h"
class Action_Principal : public Action {
  public:
    Action_Principal();
  private:
    //Lapack_Diag Principal_;
    //Lapack_General General_;
    bool doRotation_;
    AtomMask mask_;

    int init();
    int setup();
    int action();
};
#endif
