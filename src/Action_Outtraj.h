#ifndef INC_ACTION_OUTTRAJ_H
#define INC_ACTION_OUTTRAJ_H
// Outtraj
#include "Action.h"
#include "TrajectoryFile.h"

class Outtraj: public Action {
    TrajectoryFile outtraj;
    double max;
    double min;
    DataSet *Dset;
  public:
    Outtraj();
    ~Outtraj();

    int init();
    //int setup();
    int action();
    void print();
};
#endif
