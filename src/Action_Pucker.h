#ifndef INC_ACTION_PUCKER_H
#define INC_ACTION_PUCKER_H
// Class: Action_Pucker
/// Calculate the ring pucker given 5 atom masks.
#include "Action.h"
class Action_Pucker: public Action {
  public:
    Action_Pucker();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Pucker(); }
    static void Help();
  private:
    DataSet *puck_;
    AtomMask M1_;
    AtomMask M2_;
    AtomMask M3_;
    AtomMask M4_;
    AtomMask M5_;
    enum PmethodType { ALTONA=0, CREMER };
    PmethodType puckerMethod_;
    bool amplitude_;
    bool useMass_;
    bool range360_;
    double offset_;

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
};
#endif  
