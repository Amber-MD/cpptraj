#ifndef INC_ACTION_AVERAGE_H
#define INC_ACTION_AVERAGE_H
#include "Action.h"
#include "ActionFrameCounter.h"
/// Sum up all coordinates and print the averaged coords in given format.
class Action_Average: public Action, ActionFrameCounter {
  public:
    Action_Average();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Average(); }
    void Help() const;
    ~Action_Average();
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
    int SyncAction();

    int ensembleNum_;
    int debug_;
    AtomMask Mask1_;
    Frame* AvgFrame_;
    Topology AvgParm_;
    ArgList trajArgs_;
    int Natom_;
    int Nframes_;
    std::string avgfilename_;
    DataSet* crdset_;         ///< DataSet to save avg coords to.
};
#endif  
