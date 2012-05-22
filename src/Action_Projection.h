#ifndef INC_ACTION_PROJECTION_H
#define INC_ACTION_PROJECTION_H
#include "Action.h"
#include "ModesInfo.h"
/// project snapshots on normal modes
class Action_Projection : public Action {
  public:
    Action_Projection();
  private:
    int init();
    int setup();
    int action();

    ModesInfo::modesType type_;
    ModesInfo modinfo_;
    int beg_;
    int end_;
    int start_;
    int stop_;
    int offset_;
    std::vector<double> sqrtmasses_;
    AtomMask mask_;
    CpptrajFile outfile_;
};
#endif
