#ifndef INC_ACTION_AVERAGE_H
#define INC_ACTION_AVERAGE_H
/// Class: Average
/// Action that will sum up all given coordinates and print the averaged coords
/// for the given format.
#include "Action.h"
class Average: public Action {
    AtomMask Mask1;
    Frame *AvgFrame;
    AmberParm *AvgParm;
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
