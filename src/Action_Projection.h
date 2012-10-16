#ifndef INC_ACTION_PROJECTION_H
#define INC_ACTION_PROJECTION_H
#include "Action.h"
#include "DataSet_Modes.h"
/// project snapshots on normal modes
class Action_Projection : public Action {
  public:
    Action_Projection();
  private:
    int init();
    int setup();
    int action();

    typedef std::vector<DataSet*> Darray;

    Darray project_;
    DataSet_Modes modinfo_;
    int beg_;
    int end_;
    int start_;
    int stop_;
    int offset_;
    std::vector<double> sqrtmasses_;
    AtomMask mask_;
};
#endif
