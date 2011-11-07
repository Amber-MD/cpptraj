#ifndef INC_ACTION_PTRAJACTION_H
#define INC_ACTION_PTRAJACTION_H
// Class: PtrajAction
/// Wrapper for ptraj functions in ptraj_actions.c
#include "Action.h"
#include "ptraj_actionheader.h"
class PtrajAction: public Action {
    actionInformation *actioninfo;
    argStackType *argumentStack;
    double *x_coord;
    double *y_coord;
    double *z_coord;
  public:
    PtrajAction();
    ~PtrajAction();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
