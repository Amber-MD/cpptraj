#ifndef INC_ACTION_RUNNINGAVG_H
#define INC_ACTION_RUNNINGAVG_H
#include "Action.h"
// Class: RunningAvg
/// Replace current frame with running average over N frames. 
class RunningAvg: public Action {
    int Nwindow;               ///< Size of the running average
    int frameThreshold;        ///< Frame above which averaging should start, Nwindow-1
    int currentWindow;         ///< Current Position in FrameCoords
    std::vector<Frame> Window; ///< Hold coords for Nwindow frames
    int windowNatom;           ///< # of atoms in each window
    Frame avgFrame;            ///< Frame to hold averaged coords.

  public:
    RunningAvg();
    ~RunningAvg();

    void SeparateInit(int, int);

    int init();
    int setup();
    int action();
    //void print();
};
#endif  
