#ifndef INC_ACTION_LESSPLIT_H
#define INC_ACTION_LESSPLIT_H
#include "Action.h"
#include "Trajout_Single.h"
/// Split LES frame/top into normal frames/tops.
class Action_LESsplit : public Action {
  public:
    Action_LESsplit() : lesAverage_(false), lesSplit_(false), lesParm_(0) {}
    ~Action_LESsplit();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_LESsplit(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
#   ifdef MPI
    int ParallelActionInit(Parallel::Comm const&);
#   endif
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    bool lesAverage_;           ///< If true, create LES average
    bool lesSplit_;             ///< If true, split LES frames
    typedef std::vector<AtomMask> MaskArray;
    MaskArray lesMasks_;        ///< Mask of each individual LES frame
    typedef std::vector<Trajout_Single*> Tarray;
    Tarray lesTraj_;            ///< Output trajectories for split LES
    std::string splitfilename_; ///< Split LES output traj prefix
    Trajout_Single avgTraj_;    ///< Output trajectory for LES average.
    ArgList trajArgs_;          ///< Split LES output trajectory arguments.
    Frame lesFrame_;            ///< Frame for LES split
    Frame avgFrame_;            ///< Frame for LES average
    Topology* lesParm_;         ///< Topology for LES split/average
#   ifdef MPI
    Parallel::Comm trajComm_;
#   endif
};
#endif
