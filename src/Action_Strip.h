#ifndef INC_ACTION_STRIP_H
#define INC_ACTION_STRIP_H
#include "Action.h"
// Class: Action_Strip
/// Used to remove atoms from the state.
class Action_Strip: public Action {
  public:
    Action_Strip();
    ~Action_Strip();
    void print() {}
  private:
    int init();
    int setup();
    int action();

    Topology *oldParm_;
    Topology *newParm_;
    Frame newFrame_;
    std::string prefix_;
    AtomMask M1_;
    bool removeBoxInfo_;
};
// Class: Action_Unstrip
/// Signals to ActionList that the original traj parm should be restored.
class Action_Unstrip: public Action {
  public:
    Action_Unstrip() {}
  private:
    int init()   {return 0;}
    int setup()  {return 2;}
    int action() {return 2;}
    void print() {}
};
#endif  
