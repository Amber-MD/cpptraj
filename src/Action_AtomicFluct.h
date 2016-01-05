#ifndef INC_ACTION_ATOMICFLUCT_H
#define INC_ACTION_ATOMICFLUCT_H
#include "Action.h"
#include "ActionFrameCounter.h"
class Action_AtomicFluct : public Action, ActionFrameCounter {
  public :
    Action_AtomicFluct();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_AtomicFluct(); }
    void Help() const;
  private :
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    int SyncAction();
    void Print();

    enum outputType { BYATOM = 0, BYRES, BYMASK };

    Frame SumCoords_;         ///< Hold the average coordinates.
    Frame SumCoords2_;        ///< Hold the variance of coordinates.
    Frame Cross_;             ///< Hold cross-terms for calculating covariance.
    CharMask Mask_;
    int sets_;
    bool bfactor_;
    bool calc_adp_;
    CpptrajFile* adpoutfile_;
    std::string outfilename_;
    Topology *fluctParm_;
    outputType outtype_;
    DataSet* dataout_;
};
#endif
