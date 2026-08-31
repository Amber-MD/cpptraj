#ifndef INC_ACTION_SURFACETENSION_H
#define INC_ACTION_SURFACETENSION_H
#include <vector>
#include "Action.h"
#include "AtomMask.h"
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
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif

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
};
#endif
