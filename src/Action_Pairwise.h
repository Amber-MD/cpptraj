#ifndef INC_ACTIONS_PAIRWISE_H
#define INC_ACTIONS_PAIRWISE_H
/// Class: Pairwise 
/// Action to compare pairs of atoms. 
#include "Action.h"
class Pairwise: public Action {
    AtomMask Mask0;
    bool *skipv;
    int *natexidx;
    bool hasExclusion;
  public:
    Pairwise();
    ~Pairwise();

    int init();
    int setup();
    int action();
};
#endif  
