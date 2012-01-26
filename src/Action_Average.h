#ifndef INC_ACTION_AVERAGE_H
#define INC_ACTION_AVERAGE_H
#include "Action.h"
// Class: Average
/// Sum up all coordinates and print the averaged coords in given format.
class Average: public Action {
    AtomMask Mask1;
    Frame *AvgFrame;
    AmberParm *AvgParm;
    ArgList trajArgs;
    bool parmStripped;
    int Natom;
    int Nframes;
    int start;
    int stop;
    int offset;
    int targetFrame;
    char *avgfilename;

    int init();
    int setup();
    int action();
    void print();
  public:
    Average();
    ~Average();
};
#endif  
