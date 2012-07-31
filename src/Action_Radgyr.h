#ifndef INC_ACTION_RADGYR_H
#define INC_ACTION_RADGYR_H
#include "Action.h"
// Class: Action_Radgyr
/// Action to calculate the radius of gyration of atoms within a mask.
class Action_Radgyr: public Action {
  public:
    Action_Radgyr();
  private:
    int init();
    int setup();
    int action();

    DataSet* rog_;
    DataSet* rogmax_;
    AtomMask Mask1_;
    bool calcRogmax_;
};
#endif
