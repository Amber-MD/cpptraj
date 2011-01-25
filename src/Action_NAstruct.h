#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// NAstruct
#include <vector>
#include "Action.h"
#include "AxisType.h"

class NAstruct: public Action {
    // Variables
    std::vector<AxisType*> RefCoords;
    std::vector<AxisType*> BaseAxes;
    std::vector<AxisType*> BasePairAxes;
    std::vector<AxisType*> ExpFrames;
    std::vector<AtomMask*> ExpMasks;
    std::vector<int> BasePair;
    int Nbp;
    int Nbp2;
    int Nbases;
    double HBcut2;
    double Ocut2;
    AxisType ExpFrame;
    AxisType RefFrame;
    // Init Args
    std::list<int> *resRange;
    char *outFilename;
    // Functions
    void ClearLists();
    bool GCpair(AxisType *, AxisType *);
    bool ATpair(AxisType *, AxisType *);
    bool basesArePaired(AxisType *, AxisType *);
    int determineBasePairing();
    int setupBasePairAxes();
    int setupBaseAxes(Frame *);
  public:
    NAstruct();
    ~NAstruct();

    int init();
    int setup();
    int action();
};
#endif
