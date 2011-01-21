#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// NAstruct
#include "Action.h"
#include <vector>

class NAstruct: public Action {
    enum NAbaseType { UNKNOWN_BASE, DA, DT, DG, DC, RA, RC, RG, RU };
    struct AxisType {
      Frame *F;
      AmberParm::NAME *Name;
      // Use ucell[9] and recip[6] in Frame for rot mat./trans vectors
      //double R[9];
      //double trans[6];
    };
    std::vector<AxisType*> RefCoords;
    std::vector<AxisType*> BaseAxes;
    std::vector<Frame*> ExpFrames;
    std::vector<AtomMask*> ExpMasks;
    typedef int BPTYPE[2];
    BPTYPE *BasePair;
    int Nbp;
    int Nbases;
    Frame *REF_TEMP;
    Frame *EXP_TEMP;

    std::list<int> *resRange;
    char *outFilename;

    AxisType *AllocAxis(int);
    void FreeAxis( AxisType * );
    void ClearLists();
    NAbaseType ID_base(char*);
    AxisType *getRefCoords( NAbaseType);
    AxisType *principalAxes();
    int determineBasePairing();
    // DEBUG
    void AxisToPDB(PtrajFile *, AxisType *, int, int *);
  public:
    NAstruct();
    ~NAstruct();

    int init();
    int setup();
    int action();
};
#endif
