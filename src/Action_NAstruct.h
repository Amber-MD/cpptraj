#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// NAstruct
#include "Action.h"
#include <vector>

class NAstruct: public Action {
    enum NAbaseType { UNKNOWN_BASE, DA, DT, DG, DC, RA, RC, RG, RU };
    static const char NAbaseName[][5];
    struct AxisType {
      Frame *F;
      AmberParm::NAME *Name;
      NAbaseType ID;
      //double R[9];
      //double trans[6];
    };
    std::vector<AxisType*> RefCoords;
    std::vector<AxisType*> BaseAxes;
    std::vector<AxisType*> BasePairAxes;
    std::vector<Frame*> ExpFrames;
    std::vector<AtomMask*> ExpMasks;
    typedef int BPTYPE[2];

    BPTYPE *BasePair;
    int Nbp;
    int Nbases;
    double HBcut2;
    double Ocut2;
    Frame *Axes1;
    Frame *Axes2;

    std::list<int> *resRange;
    char *outFilename;

    AxisType *AllocAxis(int);
    AxisType *AxisCopy( AxisType * );
    void AxisToFrame( AxisType *, Frame * );
    void FreeAxis( AxisType * );
    void ClearLists();
    NAbaseType ID_base(char*);
    AxisType *getRefCoords( NAbaseType);
    AxisType *principalAxes();
    bool GCpair(AxisType *, AxisType *);
    bool ATpair(AxisType *, AxisType *);
    bool basesArePaired(AxisType *, AxisType *);
    int determineBasePairing();
    void flipAxesYZ( Frame * );
    int setupBasePairAxes();
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
