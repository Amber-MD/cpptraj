#ifndef INC_ANALYSIS_CALCDIFFUSION_H
#define INC_ANALYSIS_CALCDIFFUSION_H
#include "Analysis.h"
#include "DiffusionResults.h"
/// Calculate diffusion from unwrapped coordinates using multiple time origins 
class Analysis_CalcDiffusion : public Analysis {
  public:
    Analysis_CalcDiffusion();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_CalcDiffusion(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
#   ifdef MPI
    bool IsParallel() const { return true; }
#   endif
  private:
    DataSet_Coords* TgtTraj_; ///< Coordinates to calculate diffusion for
    int debug_;
    int maxlag_;              ///< Maximum number of frames to use for each origin
    double time_;             ///< Time between frames in ps
    AtomMask mask1_;          ///< Atoms to track diffusion for.
    DataSet* avg_x_;  ///< Hold average diffusion in X direction each frame
    DataSet* avg_y_;  ///< Hold average diffusion in Y direction each frame
    DataSet* avg_z_;  ///< Hold average diffusion in Z direction each frame
    DataSet* avg_r_;  ///< Hold average MSD each frame
    DataSet* avg_a_;  ///< Hold average distance each frame
    Cpptraj::DiffusionResults results_;
#   ifdef MPI
    Parallel::Comm trajComm_;
#   endif
};
#endif
