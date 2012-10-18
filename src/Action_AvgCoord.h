#ifndef INC_ACTION_AVGCOORD_H
#define INC_ACTION_AVGCOORD_H
/// Class: Action_AvgCoord
/// Action to calculate the overall average of all atomic coords. 
#include "Action.h"
class Action_AvgCoord: public Action {
  public:
    Action_AvgCoord();
    ~Action_AvgCoord();

    void print() {}
  private:
    int init();
    int setup();
    int action();

    bool calcMagnitude_;
    bool useMass_;
    AtomMask Mask_;
    CpptrajFile outfile_;
};
#endif
