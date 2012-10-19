#ifndef INC_ACTION_CENTER_H
#define INC_ACTION_CENTER_H
/// Class: Action_Center
/// Action to center coordinates to coord origin or box origin.
#include "Action.h"
class Action_Center: public Action {
  public:
    Action_Center();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Center(); }
    static void Help();

    void print() {}
  private:
    int init();
    int setup();
    int action();

    AtomMask Mask_;
    bool origin_;
    bool useMass_;
};
#endif
