#ifndef INC_ACTION_OUTTRAJ_H
#define INC_ACTION_OUTTRAJ_H
// Action_Outtraj
#include "Action.h"
#include "Trajout.h"
/// Write out a trajectory inside the ActionList
class Action_Outtraj: public Action {
  public:
    Action_Outtraj();
  private:
    int init();
    int action();
    void print();

    Trajout outtraj_;
    double max_;
    double min_;
    DataSet* Dset_;
};
#endif
