#ifndef INC_ACTION_AVERAGE_H
#define INC_ACTION_AVERAGE_H
#include "Action.h"
#include "ActionFrameCounter.h"
#include "Trajout_Single.h"
/// Sum up all coordinates and print the averaged coords in given format.
class Action_Average: public Action, ActionFrameCounter {
  public:
    Action_Average();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Average(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction(Parallel::Comm const&);
#   endif
    void Print();

    int debug_;
    AtomMask Mask1_;          ///< Mask of atoms to average.
    Frame AvgFrame_;          ///< Hold averaged coordinates.
    Topology AvgParm_;        ///< Hold topology corresponding to averaged coordinates.
    int Nframes_;             ///< Number of frames in the average.
    Trajout_Single outtraj_;  ///< File to write avg coords to.
    DataSet* crdset_;         ///< DataSet to save avg coords to.
};
#endif  
