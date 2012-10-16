#ifndef INC_ACTION_OUTTRAJ_H
#define INC_ACTION_OUTTRAJ_H
// Action_Outtraj
#include "Action.h"
#include "TrajectoryFile.h"
/// Write out a trajectory inside the ActionList
class Action_Outtraj: public Action {
  public:
    Action_Outtraj();
  private:
    int init();
    int setup() { return 0; }
    int action();
    void print();

    TrajectoryFile outtraj_;
    std::vector<double> Max_;
    std::vector<double> Min_;
    std::vector<DataSet*> Dsets_;
    DataSet* maxmin_;
};
#endif
