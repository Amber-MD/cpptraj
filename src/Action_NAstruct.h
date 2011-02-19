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
    std::vector<int> BasePair;           // Base1,Base2,IsAnti,...
    int Nbp;
    int Nbases;
    double HBcut2;
    double Ocut2;
    AxisType ExpFrame;
    AxisType RefFrame;
    // Datasets
    std::vector<DataSet*> SHEAR;
    std::vector<DataSet*> STRETCH;
    std::vector<DataSet*> STAGGER;
    std::vector<DataSet*> BUCKLE;
    std::vector<DataSet*> PROPELLER;
    std::vector<DataSet*> OPENING;
    std::vector<DataSet*> SHIFT;
    std::vector<DataSet*> SLIDE;
    std::vector<DataSet*> RISE;
    std::vector<DataSet*> TILT;
    std::vector<DataSet*> ROLL;
    std::vector<DataSet*> TWIST;
    int Nframe;                          // Keep track of # frames for print() function           
    // Init Args
    Range resRange;
    char *outFilename;
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
