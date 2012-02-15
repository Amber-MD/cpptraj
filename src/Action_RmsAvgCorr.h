#ifndef INC_ACTION_RMSAVGCORR_H
#define INC_ACTION_RMSAVGCORR_H
#include "Action.h"
#include "Action_Rmsd.h"
#include "Action_RunningAvg.h"
// Class: RmsAvgCorr
/// Calculate rmsd using running avg structures 
class RmsAvgCorr: public Action {
    Rmsd rmsdaction;
    ArgList rmsdArglist;
    DataSetList localdata;
    FrameList ReferenceFrames;
    DataSet *Ct;
    int parmNatom;
  public:
    RmsAvgCorr();
    ~RmsAvgCorr();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
