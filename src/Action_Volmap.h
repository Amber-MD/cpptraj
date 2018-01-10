#ifndef INC_ACTION_VOLMAP_H
#define INC_ACTION_VOLMAP_H
#include "Action.h"
#include "DataSet_GridFlt.h"
class Action_Volmap : public Action {
  public:
    Action_Volmap();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Volmap(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
    void RawHelp() const;

    /// grid resolutions
    double dx_, dy_, dz_;
    /// minimum values in the x-, y-, and z-dimensions
    double xmin_, ymin_, zmin_;
    int Nframes_;           ///< Number of frames we analyzed so we can average at the end
    bool setupGridOnMask_;  ///< If true, set up the grid on first frame based on centermask.
    bool spheremode_;       ///< If true, grid points farther than rhalf^2 will be skipped
    AtomMask centermask_;   ///< Mask to center the grid on
    AtomMask densitymask_;  ///< Max of atoms to grid.
    DataSet_GridFlt* grid_; ///< Hold the grid.
    DataSet* total_volume_; ///< Hold total grid volume.
    CpptrajFile* peakfile_; ///< file name with the peak locations as Carbons in XYZ file format
    double peakcut_;        ///< The value below which to ignore all peaks
    std::vector<float> halfradii_; ///< 1/2 the atomic radii of each atom in the gridded selection
    double buffer_;         ///< Clearance between the edges of our grid and centermask_
    double radscale_;       ///< The scaling factor to divide all radii by
    double stepfac_;        ///< Factor for determining how many steps to smear Gaussian
    static const double sqrt_8_pi_cubed;
#   ifdef _OPENMP
    typedef std::vector< Grid<float> > Garray;
    Garray GRID_THREAD_;
    void CombineGridThreads();
#   endif
};
#endif
