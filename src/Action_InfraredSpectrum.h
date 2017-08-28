#ifndef INC_ACTION_INFRAREDSPECTRUM_H
#define INC_ACTION_INFRAREDSPECTRUM_H
#include "Action.h"
#include "DataSet_Vector.h"
/// <Enter description of Action_InfraredSpectrum here>
class Action_InfraredSpectrum : public Action {
  public:
    Action_InfraredSpectrum();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_InfraredSpectrum(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    typedef DataSet_Vector VelArray;
    VelArray* Vel_;         ///< Hold velocity*charge for selected atoms at each frame.
    AtomMask mask_;         ///< Atoms to calculate VAC fn for.
    Topology* currentTop_;
    int maxLag_;            ///< Maximum lag to calculate VAC fn for.
    int previousNselected_; ///< Used to check if selected # atoms has changed.
    bool useFFT_;           ///< Use FFT to calculate VAC functions
};
#endif
