#ifndef INC_ACTION_RMSAVGCORR_H
#define INC_ACTION_RMSAVGCORR_H
#include "Action.h"
#include "Action_Rmsd.h"
#include "Action_RunningAvg.h"
// Class: RmsAvgCorr
/// Calculate rmsd using running avg structures 
class RmsAvgCorr: public Action {
    char *rmsmask;
    FrameList ReferenceFrames;
    DataSet *Ct;
    int parmNatom;
    int maxframes;
  public:
    RmsAvgCorr();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
