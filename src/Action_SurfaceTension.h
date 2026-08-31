#ifndef INC_ACTION_SURFACETENSION_H
#define INC_ACTION_SURFACETENSION_H
#include <vector>
#include "Action.h"
#include "AtomMask.h"
class DataSet_Mesh;
/// Capillary-wave surface tension of a liquid slab (Cartesian normal).
/** Instantaneous upper and lower interfaces are either a Willard-Chandler
  * isosurface of a Gaussian-smoothed number-density field (default) or an
  * ITIM-style per-column min/max of <mask>, split at mid-box along the normal.
  * Height fluctuations in the interface plane are Fourier transformed with
  * the numpy convention
  *   h_q = (1 / N₁N₂) Σ h(t1,t2) exp(−i q · r)
  * Capillary-wave theory including bending rigidity κ is
  *   ⟨|h_q|²⟩ = k_B T / (A (γ q² + κ q⁴))
  * so γ (mN/m) is obtained from the small-q plateau of q² S(q), and κ (in kT)
  * from the slope of 1/(q² S) vs q² on the same window.
  * The slab normal is x, y, or z (default z); gridspacing/sigmaxy apply in the
  * interface plane and dz/sigmaz along the normal. Lateral box lengths are
  * held fixed (NVT). DataSets are always created; files are written only when
  * the matching *out keyword is given. MPI: packed ReduceMaster of |h_q|².
  * \author Nathan D Levinzon <ndlevinzon@gmail.com>
  */
class Action_SurfaceTension : public Action {
  public:
    Action_SurfaceTension();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_SurfaceTension(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
#   ifdef MPI
    /// Sum Fourier accumulators onto the master rank.
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif

    /// Allocate density / height / power arrays for the lateral grid (and nz).
    int AllocateGrid(int, int, int);
    /// Build interfaces and accumulate spectra for one frame.
    /** \return 0 OK, 1 skip frame, 2 fatal (grid or lateral box changed). */
    int ProcessFrame(Frame const&, double, double, double);
    /// Finish one complete nblock window (γ, κ, and roughness).
    int FinishBlock();

    /// Instantaneous-interface definition.
    enum IfaceType { WILLARD = 0, ITIM };
    /// Cartesian slab normal.
    enum NormalAxis { AXIS_X = 0, AXIS_Y = 1, AXIS_Z = 2 };

    AtomMask Mask_;                 ///< Atoms used for the number-density / ITIM field
    IfaceType iface_;               ///< Willard-Chandler isosurface or ITIM min/max
    NormalAxis normal_;             ///< Slab normal (default z)

    double temp_;                   ///< Temperature (K)
    double gridspacing_;            ///< Target bin spacing in the interface plane (Å)
    double dz_;                     ///< Target bin spacing along the slab normal (Å)
    double sigma_xy_;               ///< Gaussian smoothing in the interface plane (Å)
    double sigma_z_;                ///< Gaussian smoothing along the slab normal (Å)
    double bulk_halfwidth_;         ///< Half-width around slab center for ρ_bulk (Å)
    double threshold_frac_;         ///< Interface is this fraction of ρ_bulk
    double qmin_;                   ///< Fit-window minimum |q| (Å⁻¹)
    double qmax_;                   ///< Fit-window maximum |q| (Å⁻¹)
    double lx_user_;                ///< Optional fixed box Lx; < 0 means use the box
    double ly_user_;                ///< Optional fixed box Ly; < 0 means use the box
    double lz_user_;                ///< Optional fixed box Lz; < 0 means use the box
    int nblock_;                    ///< Frames per uncertainty block; 0 disables
    int debug_;                     ///< Debug level from ActionInit

    // ----- Spectrum vs q; filled in Print() --------------------------------
    DataSet_Mesh* S_;               ///< Combined S(q) = ⟨|h_q|²⟩ (Å²)
    DataSet_Mesh* S_top_;           ///< Upper-interface S(q)
    DataSet_Mesh* S_bot_;           ///< Lower-interface S(q)
    DataSet_Mesh* q2S_;             ///< Combined q² S(q)
    DataSet_Mesh* q2S_top_;         ///< Upper q² S(q)
    DataSet_Mesh* q2S_bot_;         ///< Lower q² S(q)
    DataSet_Mesh* gammaq_;          ///< Apparent γ(q) (mN/m), combined
    DataSet_Mesh* gammaq_top_;      ///< Apparent γ(q), upper
    DataSet_Mesh* gammaq_bot_;      ///< Apparent γ(q), lower
    DataSet_Mesh* kappaq_;          ///< Apparent κ(q) (kT), combined
    DataSet_Mesh* kappaq_top_;      ///< Apparent κ(q), upper
    DataSet_Mesh* kappaq_bot_;      ///< Apparent κ(q), lower
    // ----- Per-frame roughness / bulk density ------------------------------
    DataSet* wtop_;                 ///< Upper RMS roughness w (Å)
    DataSet* wbot_;                 ///< Lower RMS roughness w (Å)
    DataSet* wmean_;                ///< Mean of upper and lower w (Å)
    DataSet* rhobulk_;              ///< Bulk number density (Å⁻³)
    // ----- Per-block results (only if nblock > 0) --------------------------
    DataSet* block_gamma_;          ///< Block γ (mN/m)
    DataSet* block_kappa_;          ///< Block κ (kT)
    DataSet* block_wmean_;          ///< Block mean roughness (Å)
    DataSet* block_wtop_;           ///< Block upper roughness (Å)
    DataSet* block_wbot_;           ///< Block lower roughness (Å)

    std::vector<double> density_;   ///< n1×n2×nz number density (Å⁻³); unused for ITIM
    std::vector<double> h_upper_;   ///< n1×n2 upper height field along the normal (Å)
    std::vector<double> h_lower_;   ///< n1×n2 lower height field along the normal (Å)
    std::vector<double> total_power_;  ///< Accumulated |h_q|², both surfaces
    std::vector<double> top_power_;    ///< Accumulated |h_q|², upper
    std::vector<double> bottom_power_; ///< Accumulated |h_q|², lower
    std::vector<double> block_power_;  ///< Combined |h_q|² for the open block
    std::vector<double> q_grid_;       ///< |q| for each Fourier mode (Å⁻¹)

    int nx_;                        ///< Bins along lateral axis 1
    int ny_;                        ///< Bins along lateral axis 2
    int nz_;                        ///< Density bins along the slab normal
    double Lt1_ref_;                ///< Lateral length 1 from the first good frame (Å)
    double Lt2_ref_;                ///< Lateral length 2 from the first good frame (Å)
    bool grid_ready_;               ///< True after the first successful frame

    int n_frames_;                  ///< Frames that contributed to the spectra
    int n_surfaces_;                ///< 2 × n_frames_ (upper + lower)
    int n_skipped_;                 ///< Frames skipped (no interface / no box)
    int n_blocks_;                  ///< Completed uncertainty blocks
    int block_surface_count_;       ///< Surfaces accumulated in the open block
    int block_frame_count_;         ///< Frames accumulated in the open block
    double block_w_sum_;            ///< Running sum of mean w in the open block
    double block_wtop_sum_;         ///< Running sum of upper w
    double block_wbot_sum_;         ///< Running sum of lower w
};
#endif
