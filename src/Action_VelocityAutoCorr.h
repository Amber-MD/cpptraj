#ifndef INC_ACTION_VELOCITYAUTOCORR_H
#define INC_ACTION_VELOCITYAUTOCORR_H
#include "Action.h"
#include "DataSet_Vector.h"
class Action_VelocityAutoCorr : public Action {
  public:
    Action_VelocityAutoCorr();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_VelocityAutoCorr(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int ParallelPreviousFramesRequired() const;
    int ParallelPreloadFrames(FArray const&);
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print();

    typedef DataSet_Vector Varray;
    typedef std::vector<Varray> VelArray;
    VelArray Vel_;         ///< Hold velocity info for each selected atom at each frame.
    AtomMask mask_;        ///< Atoms to calculate VAC fn for.
    Frame previousFrame_;  ///< Hold previous frame coords (!useVelInfo only)
    CpptrajFile* diffout_; ///< File to write diffusion constants to (or STDOUT)
    DataSet* VAC_;         ///< Hold values of the velocity auto-correlation function
    DataSet* diffConst_;   ///< Hold value of diffusion constant in E-5 cm^2/s
    double tstep_;         ///< Time step between frames
    int maxLag_;           ///< Maximum lag to calculate VAC fn for.
    bool useVelInfo_;      ///< If true use actual velocities in frame if present
    bool useFFT_;          ///< Use FFT to calculate VAC functions
    bool normalize_;       ///< Normalize VAC fn to 1.0 
};
#endif
