#ifndef INC_ACTION_NASTRUCT_H
#define INC_ACTION_NASTRUCT_H
// NAstruct
#include "Action.h"

class NAstruct: public Action {
    enum NAbaseType { UNKNOWN_BASE, DA, DT, DG, DC, RA, RC, RG, RU };
    typedef char NAME[5];
    struct AxisType {
      Frame *F;
      NAME *Name;
      // Use ucell and recip in Frame
      //double R[9];
      //double trans[6];
    };
    std::list<AxisType*> RefCoords;

    std::list<int> *resRange;
    char *outFilename;

    NAbaseType ID_base(char*);
    AxisType *getRefCoords( NAbaseType);
  public:
    NAstruct();
    ~NAstruct();

    int init();
    int setup();
    int action();
};
#endif
