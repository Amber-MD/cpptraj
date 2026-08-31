#ifndef INC_ACTION_SURFACETENSION_H
#define INC_ACTION_SURFACETENSION_H
#include <vector>
#include "Action.h"
#include "AtomMask.h"
<<<<<<< HEAD
<<<<<<< HEAD
#include "PubFFT.h"
class DataSet_Mesh;
class CpptrajFile;
/// Capillary-wave surface tension of a liquid slab (Cartesian normal).
/** Instantaneous interfaces are a Willard–Chandler isosurface of a
  * Gaussian-smoothed number-density field (default) or ITIM min/max.
  * Default is a two-interface slab (vacuum or a second phase on both
  * sides of the film). nsurf 1 is a single interface. An optional
  * second mask supplies the lower surface (leaflet / liquid–liquid).
  * Lateral lengths L₁, L₂ come from the unit cell unless lx/ly/lz are
  * set. Area A = L₁ L₂.
  *
  * Discrete Fourier modes of the height (FFTW / numpy unnormalized
  * transform, then divide by N):
  *   h_q = (1/(N₁ N₂)) Σ_{n₁,n₂} (h − ⟨h⟩) exp(−i q · r)
  *   q = (2π n₁/L₁, 2π n₂/L₂),   S(q) ≡ ⟨|h_q|²⟩
  *
  * Capillary-wave theory (Helfrich):
  *   S(q) = k_B T / [A (γ q² + κ q⁴)]
  * Small-q (κ → 0):  q² S(q) → k_B T / (A γ), so
  *   γ = k_B T / (A ⟨q² S(q)⟩)     on [q_min, q_max]
  * reported in mN/m (×1000 from N/m). If q_min is omitted it is the
  * fundamental wavevector 2π / max(L₁, L₂).
  *
  * Bending modulus from the linear form
  *   1/(q² S) = (A / k_B T) (γ + κ q²)
  * intercept → γ, slope → κ / k_B T.
  *
  * \author Nathan D. Levinzon <ndlevinzon@gmail.com>
  *
  * References (for interested parties or lowly PhD students like me):
  *   Willard–Chandler instantaneous interface
  *     A. P. Willard and D. Chandler, J. Phys. Chem. B 114, 1954 (2010).
  *     https://doi.org/10.1021/jp909219k
  *   ITIM (Identification of Truly Interfacial Molecules)
  *     L. B. Pártay, G. Hantal, P. Jedlovszky, Á. Vincze, and G. Horvai,
  *     J. Comput. Chem. 29, 945 (2008).
  *     https://doi.org/10.1002/jcc.20852
  *   Capillary-wave theory
  *     F. P. Buff, R. A. Lovett, and F. H. Stillinger,
  *     Phys. Rev. Lett. 15, 621 (1965).
  *     https://doi.org/10.1103/PhysRevLett.15.621
  *   Height-fluctuation surface tension in MD
  *     S. W. Sides, G. S. Grest, and M.-D. Lacasse,
  *     Phys. Rev. E 60, 6708 (1999).
  *     https://doi.org/10.1103/PhysRevE.60.6708
  *   Helfrich bending energy
  *     W. Helfrich, Z. Naturforsch. C 28, 693 (1973).
  *     https://doi.org/10.1515/znc-1973-11-1209
=======
=======
#include "PubFFT.h"
>>>>>>> 5b484999 (Implement 2-D FFT for height fields in Action_SurfaceTension)
class DataSet_Mesh;
/// Capillary-wave surface tension of a liquid slab (Cartesian normal).
/** Instantaneous upper and lower interfaces are either a Willard-Chandler
  * isosurface of a Gaussian-smoothed number-density field (default) or an
  * ITIM-style per-column min/max of <mask>, split at mid-box along the normal.
  * Height fluctuations in the interface plane are Fourier transformed with
  * cpptraj PubFFT (row-column 1-D FFTs) using the numpy convention
  *   h_q = (1 / N₁N₂) FFT2(h − ⟨h⟩)
  * Capillary-wave theory including bending rigidity κ is
  *   ⟨|h_q|²⟩ = k_B T / (A (γ q² + κ q⁴))
  * so γ (mN/m) is obtained from the small-q plateau of q² S(q), and κ (in kT)
  * from the slope of 1/(q² S) vs q² on the same window.
  * The slab normal is x, y, or z (default z); gridspacing/sigmaxy apply in the
  * interface plane and dz/sigmaz along the normal. Lateral box lengths are
  * held fixed (NVT). DataSets are always created; files are written only when
  * the matching *out keyword is given. MPI: packed ReduceMaster of |h_q|².
  * \author Nathan D Levinzon <ndlevinzon@gmail.com>
>>>>>>> c51400c7 (Add 'surftension' command to calculate capillary-wave surface tension of a liquid slab)
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
<<<<<<< HEAD
=======
    /// Sum Fourier accumulators onto the master rank.
>>>>>>> c51400c7 (Add 'surftension' command to calculate capillary-wave surface tension of a liquid slab)
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif

<<<<<<< HEAD
    int AllocateGrid(int, int, int);
    /** \return 0 OK, 1 skip frame, 2 fatal (grid or lateral box changed). */
    int ProcessFrame(Frame const&, double, double, double);
    int FinishBlock();
    void HeightPower(std::vector<double> const&, std::vector<double>&);
    /// Willard–Chandler heights from one atom set into the given density buffer.
    /** \return 0 OK, 1 skip. */
    int WillardHeights(std::vector<double> const&, std::vector<double> const&,
                       std::vector<double> const&, int,
                       double, double, double,
                       std::vector<double>&, bool, bool, double&);
    void SkipWarn(const char*);

    enum IfaceType { WILLARD = 0, ITIM };
    enum NormalAxis { AXIS_X = 0, AXIS_Y = 1, AXIS_Z = 2 };
    enum SideType { SIDE_UPPER = 0, SIDE_LOWER };

    AtomMask Mask_;                 ///< Primary density / ITIM atoms (upper if mask2)
    AtomMask Mask2_;                ///< Optional lower-surface atoms (leaflet / 2nd liquid)
    IfaceType iface_;
    NormalAxis normal_;
    SideType side_;                 ///< Which interface when nsurf == 1
    int nsurf_;                     ///< 1 or 2 instantaneous interfaces
    bool has_mask2_;
    bool do_upper_;
    bool do_lower_;
    bool qmin_specified_;

    double temp_;
    double gridspacing_;
    double dz_;
    double sigma_xy_;
    double sigma_z_;
    double bulk_halfwidth_;
    double threshold_frac_;
    double qmin_;
    double qmax_;
    double q_fundamental_;            ///< 2π / max(L₁, L₂) from the first good frame
    double lx_user_;
    double ly_user_;
    double lz_user_;
    double dt_;                     ///< Analyzed-frame spacing (ps); < 0 unused
    double blocktime_;              ///< Block length (ps); < 0 unused
    int nblock_;
    int debug_;
    int n_skip_warn_;               ///< Skip-frame warnings already printed

    DataSet_Mesh* S_;
    DataSet_Mesh* S_top_;
    DataSet_Mesh* S_bot_;
    DataSet_Mesh* q2S_;
    DataSet_Mesh* q2S_top_;
    DataSet_Mesh* q2S_bot_;
    DataSet_Mesh* gammaq_;
    DataSet_Mesh* gammaq_top_;
    DataSet_Mesh* gammaq_bot_;
    DataSet_Mesh* kappaq_;
    DataSet_Mesh* kappaq_top_;
    DataSet_Mesh* kappaq_bot_;
    DataSet* wtop_;
    DataSet* wbot_;
    DataSet* wmean_;
    DataSet* rhobulk_;
    DataSet* block_gamma_;
    DataSet* block_kappa_;
    DataSet* block_wmean_;
    DataSet* block_wtop_;
    DataSet* block_wbot_;
    CpptrajFile* summaryFile_;      ///< Optional key/value summary (summaryout)

    std::vector<double> density_;   ///< Willard field for Mask_
    std::vector<double> density2_;  ///< Willard field for Mask2_
    std::vector<double> h_upper_;
    std::vector<double> h_lower_;
    std::vector<double> total_power_;
    std::vector<double> top_power_;
    std::vector<double> bottom_power_;
    std::vector<double> block_power_;
    std::vector<double> q_grid_;

    PubFFT fft_n1_;
    PubFFT fft_n2_;
    ComplexArray fft_grid_;
    ComplexArray fft_row_;
    ComplexArray fft_col_;

    int nx_;
    int ny_;
    int nz_;
    double Lt1_ref_;
    double Lt2_ref_;
    bool grid_ready_;

    int n_frames_;
    int n_surfaces_;
    int n_skipped_;
    int n_blocks_;
    int block_surface_count_;
    int block_frame_count_;
    double block_w_sum_;
    double block_wtop_sum_;
    double block_wbot_sum_;
=======
    /// Allocate density / height / power arrays for the lateral grid (and nz).
    int AllocateGrid(int, int, int);
    /// Build interfaces and accumulate spectra for one frame.
    /** \return 0 OK, 1 skip frame, 2 fatal (grid or lateral box changed). */
    int ProcessFrame(Frame const&, double, double, double);
    /// Finish one complete nblock window (γ, κ, and roughness).
    int FinishBlock();
    /// |h_q|² of a height field via 2-D PubFFT / (nx ny).
    void HeightPower(std::vector<double> const&, std::vector<double>&);

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

    PubFFT fft_n1_;                 ///< 1-D FFT along lateral axis 1 (nx)
    PubFFT fft_n2_;                 ///< 1-D FFT along lateral axis 2 (ny)
    ComplexArray fft_grid_;         ///< nx×ny complex grid for the 2-D FFT
    ComplexArray fft_row_;          ///< Length-ny row buffer
    ComplexArray fft_col_;          ///< Length-nx column buffer

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
>>>>>>> c51400c7 (Add 'surftension' command to calculate capillary-wave surface tension of a liquid slab)
};
#endif
