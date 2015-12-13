#ifndef INC_ACTION_RUNNINGAVG_H
#define INC_ACTION_RUNNINGAVG_H
#include "Action.h"
// Class: Action_RunningAvg
/// Replace current frame with running average over N frames. 
class Action_RunningAvg: public Action {
  public:
    Action_RunningAvg();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_RunningAvg(); }
    void Help() const;
  private:
    int Nwindow_;               ///< Size of the running average
    double d_Nwindow_;          ///< For frame division (avoids constant recasting)
    int frameThreshold_;        ///< Frame above which averaging should start, Nwindow-1
    int currentWindow_;         ///< Current Position in FrameCoords
    std::vector<Frame> Window_; ///< Hold coords for Nwindow frames
    int windowNatom_;           ///< # of atoms in each window
    Frame avgFrame_;            ///< Frame to hold sum of coords in window to be avgd.
    Frame resultFrame_;         ///< Frame to hold result of averaging coords.

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}
};
#endif  
