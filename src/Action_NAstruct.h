#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// NAstruct
#include <vector>
#include "Action.h"
#include "AxisType.h"
#include "Range.h"

class NAstruct: public Action {
    // Variables
    std::vector<AxisType*> RefCoords;
    std::vector<AxisType*> BaseAxes;
    std::vector<AxisType*> BasePairAxes;
    std::vector<AxisType*> ExpFrames;
    std::vector<AtomMask*> ExpMasks;
    std::vector<int> BasePair;           // [Base1-0,Base2-0,IsAnti-0], [Base1-1...
    int Nbp;
    int Nbases;
    double HBcut2;
    double Ocut2;
    AxisType ExpFrame;
    AxisType RefFrame;
    // Datasets
    DataSetList SHEAR;
    DataSetList STRETCH;
    DataSetList STAGGER;
    DataSetList BUCKLE;
    DataSetList PROPELLER;
    DataSetList OPENING;
    DataSetList SHIFT;
    DataSetList SLIDE;
    DataSetList RISE;
    DataSetList TILT;
    DataSetList ROLL;
    DataSetList TWIST;
    int Nframe;                          // Keep track of # frames for print() function           
    // Init Args
    Range resRange;
    char *outFilename;
    char *naoutFilename;
    bool noheader;
    // Functions
    void ClearLists();
    bool GCpair(AxisType *, AxisType *);
    bool ATpair(AxisType *, AxisType *);
    bool basesArePaired(AxisType *, AxisType *);
    int determineBasePairing();
    //int setupBasePairAxes();
    int setupBaseAxes(Frame *);
    int determineBaseParameters();
    int determineBasepairParameters();
  public:
    NAstruct();
    ~NAstruct();

    int init();
    int setup();
    int action();
    void print();
};
#endif
