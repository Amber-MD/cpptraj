#ifndef INC_ACTION_RMSAVGCORR_H
#define INC_ACTION_RMSAVGCORR_H
#include "Action.h"
#include "CoordList.h"
// Class: RmsAvgCorr
/// Calculate rmsd using running avg structures 
class RmsAvgCorr: public Action {
    char *separateName;
    AtomMask Mask0;
    CoordList ReferenceCoords;
    AmberParm *ReferenceParm;
    DataSet *Ct;
    int parmNatom;
    int maxwindow;
  public:
    RmsAvgCorr();

    int init();
    int setup();
    int action();
    void print();
};
#endif  
