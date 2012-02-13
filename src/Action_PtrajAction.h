#ifndef INC_ACTION_PTRAJACTION_H
#define INC_ACTION_PTRAJACTION_H
#include "Action.h"
#include "ptraj_actions.h"
#include "ptraj_arg.h"
// Class: PtrajAction
/// Wrapper for ptraj functions in ptraj_actions.c
class PtrajAction: public Action {
    actionInformation *actioninfo;
    argStackType *argumentStack;
    double *x_coord;
    double *y_coord;
    double *z_coord;
    double ptraj_box[6];
    bool CalledSetup;
    bool coordinate_update;
  public:
    PtrajAction();
    ~PtrajAction();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
