#ifndef INC_ACTION_OUTTRAJ_H
#define INC_ACTION_OUTTRAJ_H
// Outtraj
#include "Action.h"
#include "TrajoutList.h"

class Outtraj: public Action {
    TrajoutList outtraj;
  public:
    Outtraj();
    ~Outtraj();

    int init();
    //int setup();
    int action();
};
#endif
