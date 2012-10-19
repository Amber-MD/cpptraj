#ifndef INC_ACTION_PRINCIPAL_H
#define INC_ACTION_PRINCIPAL_H
#include "Action.h"
class Action_Principal : public Action {
  public:
    Action_Principal();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Principal(); }
    static void Help();

    void print() {}
  private:
    bool doRotation_;
    bool useMass_;
    AtomMask mask_;

    int init();
    int setup();
    int action();
};
#endif
