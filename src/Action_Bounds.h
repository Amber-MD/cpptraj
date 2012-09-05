#ifndef INC_ACTION_BOUNDS_H
#define INC_ACTION_BOUNDS_H
#include "Action.h"
/// Report the min/max XYZ values for atoms in mask.
class Action_Bounds : public Action {
  public:
    Action_Bounds();
    void print();
  private:
    int init();
    int setup();
    int action();
    AtomMask mask_;
    std::string outfilename_;
    double max_[3];
    double min_[3];
};
#endif
