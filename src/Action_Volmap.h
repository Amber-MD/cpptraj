#ifndef INC_ACTION_VOLMAP_H
#define INC_ACTION_VOLMAP_H
#include "Action.h"
#include "Grid.h"
#include "SplineFxnTable.h"
#ifdef VOLMAP_DOUBLE
# define VOLMAP_DS_T DataSet_GridDbl
# define VOLMAP_T double
#else
# define VOLMAP_DS_T DataSet_GridFlt
# define VOLMAP_T float
#endif
class VOLMAP_DS_T;
/// Calculate atomic volumetric density maps from trajectory data.
/** By default the grid type used is single-precision, mostly to save space.
  * A double-precision grid can be used by compiling with the
  * VOLMAP_DOUBLE define.
  * Also by default the exp() function will be approximated with cubic spline
  * interpolation. To use the system exp() function, compile with the
  * VOLMAP_USEEXP define.
  */
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

    /// Radii to use
    enum RadiiType { UNSPECIFIED = 0, VDW, ELEMENT };
    RadiiType radiiType_;
    /// grid resolutions
    double dx_, dy_, dz_;
    /// minimum values in the x-, y-, and z-dimensions
    double xmin_, ymin_, zmin_;
    int debug_;
    int Nframes_;           ///< Number of frames we analyzed so we can average at the end
    bool setupGridOnMask_;  ///< If true, set up the grid on first frame based on centermask.
    bool spheremode_;       ///< If true, grid points farther than rhalf^2 will be skipped
    AtomMask centermask_;   ///< Mask to center the grid on
    AtomMask densitymask_;  ///< Max of atoms to grid.
    VOLMAP_DS_T* grid_;     ///< Hold the grid.
    DataSet* total_volume_; ///< Hold total grid volume.
    bool calcpeaks_;        ///< If true, calculate peaks
    DataFile* peakfile_;    ///< Optional file to write peak data to as Carbons in XYZ file format
    DataSet* peakdata_;     ///< Data set holding peak locations along with densities
    double peakcut_;        ///< The value below which to ignore all peaks
    std::vector<int> Atoms_; ///< Atoms with radii > 0.0
    std::vector<double> halfradii_; ///< 1/2 the atomic radii of each atom in the gridded selection
    double buffer_;         ///< Clearance between the edges of our grid and centermask_
    double radscale_;       ///< The scaling factor to divide all radii by
    double stepfac_;        ///< Factor for determining how many steps to smear Gaussian
    static const double sqrt_8_pi_cubed_;
    SplineFxnTable table_;
    double splineDx_;
#   ifdef _OPENMP
    typedef std::vector< Grid<VOLMAP_T> > Garray;
    Garray GRID_THREAD_;
    void CombineGridThreads();
#   endif
};
#undef VOLMAP_DS_T
#undef VOLMAP_T
#endif
