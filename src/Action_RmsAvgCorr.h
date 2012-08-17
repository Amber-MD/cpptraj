#ifndef INC_ACTION_RMSAVGCORR_H
#define INC_ACTION_RMSAVGCORR_H
#include "Action.h"
#include "CoordList.h"
// Class: Action_RmsAvgCorr
/// Calculate rmsd using running avg structures 
class Action_RmsAvgCorr: public Action {
  public:
    Action_RmsAvgCorr();

    void print();
  private:
    int init();
    int setup();
    int action();

    std::string separateName_;
    AtomMask Mask0_;
    CoordList ReferenceCoords_;
    Topology* ReferenceParm_;
    DataSet* Ct_;
    int parmNatom_;
    int maxwindow_;
};
#endif  
