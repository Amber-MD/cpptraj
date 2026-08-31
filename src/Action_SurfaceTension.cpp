// Action_SurfaceTension
// Capillary-wave surface tension of a liquid slab (Cartesian normal x, y, or z).
// Formulae, units, and references: Action_SurfaceTension.h
//
// Per-frame pipeline (after Init / Setup):
//   1. Permute Cartesian (x,y,z) into (t1, t2, n) for the chosen normal.
//   2. Wrap laterals into [0, L_α). Recenter the normal so the circular
//      mean of the film sits at L_n/2 (one shared shift if mask2 is set).
//   3. Instantaneous height h(t1, t2): Willard–Chandler isosurface or ITIM.
//   4. Roughness w = √⟨(h − ⟨h⟩)²⟩ and power |h_q|² from the DFT
//        h_q = (1/(N₁ N₂)) Σ (h − ⟨h⟩) exp(−i q · r)
//   5. Print() forms S(q) ≡ ⟨|h_q|²⟩, shell-averages, then
//        S(q) = k_B T / [A (γ q² + κ q⁴)]
//
// NVT: first good frame freezes N₁, N₂, L₁, L₂ (and N_n for Willard).
// \author Nathan D. Levinzon <ndlevinzon@gmail.com>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sys/stat.h>
#include <vector>
#include "Action_SurfaceTension.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "DataFile.h"
#include "DataSet_1D.h"
#include "DataSet_Mesh.h"
#include "CpptrajFile.h"
#include "FileName.h"
#include "Frame.h"
#ifdef _OPENMP
# include <omp.h>
#endif

// File-local helpers. Names are prefixed ST_ so they do not collide with
// other Action translation units.

/// k_B (J/K), CODATA 2018. γ uses SI: k_B T / (A_m² ⟨q² S⟩), then ×1000 → mN/m.
static const double ST_KB = 1.380649e-23;
/// Å² → m² so A in the CWT formula is consistent with γ in N/m.
static const double ST_ANG2_TO_M2 = 1.0e-20;
/// scipy.ndimage.gaussian_filter default: kernel radius = truncate × σ (pixels).
static const double ST_GAUSS_TRUNCATE = 4.0;

/// \return true if x is finite (not NaN or ±Inf). Used to reject failed crossings.
static inline bool ST_Finite(double x) {
  return (x == x) &&
         (x <  std::numeric_limits<double>::infinity()) &&
         (x > -std::numeric_limits<double>::infinity());
}

static inline double ST_NaN() {
  return std::numeric_limits<double>::quiet_NaN();
}

/// Linear index of ρ(i₁, i₂, i_n) with the normal the fastest dimension.
static inline size_t ST_Idx3(int ix, int iy, int iz, int ny, int nz) {
  return ((size_t)ix * (size_t)ny + (size_t)iy) * (size_t)nz + (size_t)iz;
}

/// Linear index of an N₁×N₂ field, row-major in the first lateral index.
static inline size_t ST_Idx2(int ix, int iy, int ny) {
  return (size_t)ix * (size_t)ny + (size_t)iy;
}

/// Round a non-negative value to the nearest integer (ties away from −∞).
static inline int ST_IRound(double x) {
  return (int)floor(x + 0.5);
}

/** Minimum-image wrap of a Cartesian coordinate into [0, L).
  * fmod of a negative argument is negative on IEEE-754; that case is lifted
  * back into the interval. x == L maps to 0 so the last bin is not empty.
  */
static double ST_Wrap(double x, double L) {
  if (L <= 0.0) return x;
  double y = fmod(x, L);
  if (y < 0.0) y += L;
  if (y >= L) y = 0.0;
  return y;
}

/** Circular mean of the film along the periodic normal, in Å.
  * Map n → θ = 2π n / L_n and take atan2(⟨sin θ⟩, ⟨cos θ⟩). A slab that
  * straddles n = 0 is then a single cluster on the circle. n2 may be 0;
  * when both masks are given they share this one mean so leaflets are not
  * pulled onto the same mid-box plane independently.
  * \return the lab-frame coordinate of that mean, in [0, L_n).
  */
static double ST_CircularSlabCenter(std::vector<double> const& n1, int n1n,
                                    std::vector<double> const* n2, int n2n,
                                    double L)
{
  if (L <= 0.0) return 0.0;
  double mean_sin = 0.0;
  double mean_cos = 0.0;
  int ntot = 0;
  for (int i = 0; i < n1n; i++) {
    double theta = Constants::TWOPI * ST_Wrap(n1[i], L) / L;
    mean_sin += sin(theta);
    mean_cos += cos(theta);
    ntot++;
  }
  if (n2 != 0) {
    for (int i = 0; i < n2n; i++) {
      double theta = Constants::TWOPI * ST_Wrap((*n2)[i], L) / L;
      mean_sin += sin(theta);
      mean_cos += cos(theta);
      ntot++;
    }
  }
  if (ntot < 1) return 0.0;
  mean_sin /= (double)ntot;
  mean_cos /= (double)ntot;
  double angle = atan2(mean_sin, mean_cos);
  if (angle < 0.0) angle += Constants::TWOPI;
  return L * angle / Constants::TWOPI;
}

/** Shift n so slab_center maps to L_n/2, then wrap into [0, L_n). */
static void ST_ApplyCircularRecenter(std::vector<double>& n, double L,
                                     double slab_center)
{
  if (n.empty() || L <= 0.0) return;
  for (size_t i = 0; i < n.size(); i++)
    n[i] = ST_Wrap(n[i] - slab_center + 0.5 * L, L);
}

/** Recenter one atom set so its circular mean along the normal lies at L_n/2. */
static void ST_CircularRecenter(std::vector<double>& n, double L) {
  if (n.empty() || L <= 0.0) return;
  ST_ApplyCircularRecenter(n, L, ST_CircularSlabCenter(n, (int)n.size(), 0, 0, L));
}

/// ASCII name of a Cartesian axis for Help / mprintf (terminals are not UTF-8).
static const char* ST_AxisName(int axis) {
  if (axis == 0) return "x";
  if (axis == 1) return "y";
  return "z";
}

/// Names of the two lateral axes for Cartesian normal ∈ {0,1,2} = {x,y,z}.
static void ST_PlaneAxes(int normal, const char*& t1, const char*& t2) {
  if (normal == 0) { t1 = "y"; t2 = "z"; }
  else if (normal == 1) { t1 = "x"; t2 = "z"; }
  else { t1 = "x"; t2 = "y"; }
}

/** Permute (L_x, L_y, L_z) into (L₁, L₂, L_n).
  *   n = z:  L₁ = L_x, L₂ = L_y, L_n = L_z
  *   n = y:  L₁ = L_x, L₂ = L_z, L_n = L_y
  *   n = x:  L₁ = L_y, L₂ = L_z, L_n = L_x
  */
static void ST_SplitBox(int normal, double Lx, double Ly, double Lz,
                        double& Lt1, double& Lt2, double& Ln)
{
  if (normal == 0) { Lt1 = Ly; Lt2 = Lz; Ln = Lx; }
  else if (normal == 1) { Lt1 = Lx; Lt2 = Lz; Ln = Ly; }
  else { Lt1 = Lx; Lt2 = Ly; Ln = Lz; }
}

/// Same permutation as ST_SplitBox, applied to one Cartesian point.
static void ST_SplitXYZ(int normal, double x, double y, double z,
                        double& t1, double& t2, double& n)
{
  if (normal == 0) { t1 = y; t2 = z; n = x; }
  else if (normal == 1) { t1 = x; t2 = z; n = y; }
  else { t1 = x; t2 = y; n = z; }
}

/** Normalized 1-D Gaussian matching scipy.ndimage._gaussian_kernel1d.
  * σ is in pixels (Å / bin width). Radius = round(4 σ). Weights are
  *   K_i = exp(−i² / (2 σ²)) / Σ_j K_j ,   i ∈ [−R, R]
  * An empty kernel (σ ≤ 0) means that axis is left unfiltered.
  */
static void ST_GaussianKernel(double sigma, std::vector<double>& kernel) {
  kernel.clear();
  if (sigma <= 0.0) return;
  int radius = ST_IRound(ST_GAUSS_TRUNCATE * sigma);
  if (radius < 0) radius = 0;
  kernel.assign((size_t)(2 * radius + 1), 0.0);
  double sigma2 = sigma * sigma;
  double sum = 0.0;
  for (int i = -radius; i <= radius; i++) {
    double v = exp(-0.5 * ((double)i * (double)i) / sigma2);
    kernel[i + radius] = v;
    sum += v;
  }
  if (sum > 0.0) {
    for (size_t i = 0; i < kernel.size(); i++)
      kernel[i] /= sum;
  }
}

/** Periodic 1-D convolution, out[i] = Σ_k K_k in[(i+k) mod N].
  * The Gaussian is even, so convolution and correlation are identical.
  */
static void ST_ConvolveWrap(std::vector<double> const& in, std::vector<double> const& kernel,
                            std::vector<double>& out)
{
  int n = (int)in.size();
  int ksize = (int)kernel.size();
  int radius = (ksize - 1) / 2;
  out.assign((size_t)n, 0.0);
  if (n == 0) return;
  if (kernel.empty() || ksize < 1) {
    out = in;
    return;
  }
  for (int i = 0; i < n; i++) {
    double acc = 0.0;
    for (int k = -radius; k <= radius; k++) {
      int j = i + k;
      j %= n;
      if (j < 0) j += n;
      acc += in[j] * kernel[k + radius];
    }
    out[i] = acc;
  }
}

/** Separable periodic Gaussian on ρ(i₁, i₂, i_n), Willard–Chandler coarse-graining.
  * σ_α in pixels (Å / Δ_α), same as scipy.ndimage.gaussian_filter. Three
  * 1-D passes (n, then t2, then t1). OpenMP: each line in a pass is
  * independent; one parallel region reuses thread-local buffers.
  */
static void ST_GaussianFilter3D(std::vector<double>& rho, int nx, int ny, int nz,
                                double sigma_x, double sigma_y, double sigma_z)
{
  std::vector<double> kx, ky, kz;
  ST_GaussianKernel(sigma_x, kx);
  ST_GaussianKernel(sigma_y, ky);
  ST_GaussianKernel(sigma_z, kz);
  if (kx.empty() && ky.empty() && kz.empty()) return;

# ifdef _OPENMP
# pragma omp parallel
# endif
  {
    std::vector<double> line, conv;
    int col;

    if (!kz.empty()) {
      line.resize((size_t)nz);
      int nxy = nx * ny;
#     ifdef _OPENMP
#     pragma omp for schedule(dynamic)
#     endif
      for (col = 0; col < nxy; col++) {
        int ix = col / ny;
        int iy = col - ix * ny;
        for (int iz = 0; iz < nz; iz++)
          line[iz] = rho[ST_Idx3(ix, iy, iz, ny, nz)];
        ST_ConvolveWrap(line, kz, conv);
        for (int iz = 0; iz < nz; iz++)
          rho[ST_Idx3(ix, iy, iz, ny, nz)] = conv[iz];
      }
    }

    if (!ky.empty()) {
      line.resize((size_t)ny);
      int nxz = nx * nz;
#     ifdef _OPENMP
#     pragma omp for schedule(dynamic)
#     endif
      for (col = 0; col < nxz; col++) {
        int ix = col / nz;
        int iz = col - ix * nz;
        for (int iy = 0; iy < ny; iy++)
          line[iy] = rho[ST_Idx3(ix, iy, iz, ny, nz)];
        ST_ConvolveWrap(line, ky, conv);
        for (int iy = 0; iy < ny; iy++)
          rho[ST_Idx3(ix, iy, iz, ny, nz)] = conv[iy];
      }
    }

    if (!kx.empty()) {
      line.resize((size_t)nx);
      int nyz = ny * nz;
#     ifdef _OPENMP
#     pragma omp for schedule(dynamic)
#     endif
      for (col = 0; col < nyz; col++) {
        int iy = col / nz;
        int iz = col - iy * nz;
        for (int ix = 0; ix < nx; ix++)
          line[ix] = rho[ST_Idx3(ix, iy, iz, ny, nz)];
        ST_ConvolveWrap(line, kx, conv);
        for (int ix = 0; ix < nx; ix++)
          rho[ST_Idx3(ix, iy, iz, ny, nz)] = conv[ix];
      }
    }
  }
}

/** Willard–Chandler crossing along one lateral column.
  * Walk from the slab mid-plane toward +n (upper) or −n (lower). The
  * interface is the first point where ρ drops through θ ρ_bulk, with θ
  * the threshold fraction. Linear interpolation between the two bins:
  *   n* = n_k + (θ ρ_bulk − ρ_k) (n_{k±1} − n_k) / (ρ_{k±1} − ρ_k)
  * \return NaN if no crossing exists (caller skips the frame).
  */
static double ST_FindCrossing(std::vector<double> const& z_grid,
                              std::vector<double> const& density,
                              double threshold, int center_index, bool upper)
{
  int nz = (int)z_grid.size();
  if (nz < 2 || center_index < 0 || center_index >= nz)
    return ST_NaN();
  if (upper) {
    for (int k = center_index; k < nz - 1; k++) {
      double rho0 = density[k];
      double rho1 = density[k + 1];
      if (rho0 >= threshold && rho1 < threshold) {
        if (fabs(rho1 - rho0) < Constants::SMALL)
          return 0.5 * (z_grid[k] + z_grid[k + 1]);
        return z_grid[k] + (threshold - rho0) * (z_grid[k + 1] - z_grid[k]) / (rho1 - rho0);
      }
    }
  } else {
    for (int k = center_index; k > 0; k--) {
      double rho0 = density[k];
      double rho1 = density[k - 1];
      if (rho0 >= threshold && rho1 < threshold) {
        if (fabs(rho1 - rho0) < Constants::SMALL)
          return 0.5 * (z_grid[k] + z_grid[k - 1]);
        return z_grid[k] + (threshold - rho0) * (z_grid[k - 1] - z_grid[k]) / (rho1 - rho0);
      }
    }
  }
  return ST_NaN();
}

/** ITIM in the probe → 0 limit, after the film has been recentered at L_n/2.
  * Each lateral column is split at mid-box:
  *   h_upper = max { n_i | n_i ≥ L_n/2 }     (need_u)
  *   h_lower = min { n_i | n_i <  L_n/2 }     (need_l)
  * An empty required half-column returns false (skip the frame). This is
  * not a finite-radius ITIM probe; it is the min/max of the mask.
  */
static bool ST_ItimMinMax(std::vector<double> const& t1,
                          std::vector<double> const& t2,
                          std::vector<double> const& n,
                          int natom, double Lt1, double Lt2, double Ln,
                          int nx, int ny,
                          bool need_u, bool need_l,
                          std::vector<double>& h_upper,
                          std::vector<double>& h_lower)
{
  double dx = Lt1 / (double)nx;
  double dy = Lt2 / (double)ny;
  double mid = 0.5 * Ln;
  size_t n2 = (size_t)nx * (size_t)ny;
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> zmax(n2, -inf), zmin(n2, inf);
  std::vector<char> has_u(n2, 0), has_l(n2, 0);
  for (int i = 0; i < natom; i++) {
    int ix = (int)floor(t1[i] / dx);
    int iy = (int)floor(t2[i] / dy);
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (ix >= nx) ix = nx - 1;
    if (iy >= ny) iy = ny - 1;
    size_t idx = ST_Idx2(ix, iy, ny);
    if (n[i] >= mid) {
      if (!has_u[idx] || n[i] > zmax[idx]) zmax[idx] = n[i];
      has_u[idx] = 1;
    } else {
      if (!has_l[idx] || n[i] < zmin[idx]) zmin[idx] = n[i];
      has_l[idx] = 1;
    }
  }
  for (size_t i = 0; i < n2; i++) {
    if (need_u && !has_u[i]) return false;
    if (need_l && !has_l[i]) return false;
    if (need_u) h_upper[i] = zmax[i];
    if (need_l) h_lower[i] = zmin[i];
  }
  return true;
}

/** Two-mask ITIM (leaflet / liquid–liquid): no mid-box split.
  *   h_upper(i₁,i₂) = max n of mask 1 in that column
  *   h_lower(i₁,i₂) = min n of mask 2 in that column
  * \return false if any required column is empty.
  */
static bool ST_ItimTwoMasks(std::vector<double> const& t1u,
                            std::vector<double> const& t2u,
                            std::vector<double> const& nu,
                            int natomu,
                            std::vector<double> const& t1l,
                            std::vector<double> const& t2l,
                            std::vector<double> const& nl,
                            int natoml,
                            double Lt1, double Lt2,
                            int nx, int ny,
                            std::vector<double>& h_upper,
                            std::vector<double>& h_lower)
{
  double dx = Lt1 / (double)nx;
  double dy = Lt2 / (double)ny;
  size_t n2 = (size_t)nx * (size_t)ny;
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> zmax(n2, -inf), zmin(n2, inf);
  std::vector<char> has_u(n2, 0), has_l(n2, 0);
  for (int i = 0; i < natomu; i++) {
    int ix = (int)floor(t1u[i] / dx);
    int iy = (int)floor(t2u[i] / dy);
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (ix >= nx) ix = nx - 1;
    if (iy >= ny) iy = ny - 1;
    size_t idx = ST_Idx2(ix, iy, ny);
    if (!has_u[idx] || nu[i] > zmax[idx]) zmax[idx] = nu[i];
    has_u[idx] = 1;
  }
  for (int i = 0; i < natoml; i++) {
    int ix = (int)floor(t1l[i] / dx);
    int iy = (int)floor(t2l[i] / dy);
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (ix >= nx) ix = nx - 1;
    if (iy >= ny) iy = ny - 1;
    size_t idx = ST_Idx2(ix, iy, ny);
    if (!has_l[idx] || nl[i] < zmin[idx]) zmin[idx] = nl[i];
    has_l[idx] = 1;
  }
  for (size_t i = 0; i < n2; i++) {
    if (!has_u[i] || !has_l[i]) return false;
    h_upper[i] = zmax[i];
    h_lower[i] = zmin[i];
  }
  return true;
}

/** Extract (t1, t2, n) for one mask. Laterals are wrapped into [0, L_α);
  * the normal is left unshifted. Recenter is applied afterwards so two
  * masks can share one circular mean.
  */
static void ST_SplitAtomCoords(Frame const& frm, AtomMask const& mask, int nax,
                               double Lx, double Ly, double Lz,
                               std::vector<double>& t1, std::vector<double>& t2,
                               std::vector<double>& n)
{
  double Lt1, Lt2, Ln;
  ST_SplitBox(nax, Lx, Ly, Lz, Lt1, Lt2, Ln);
  int natom = mask.Nselected();
  t1.resize((size_t)natom);
  t2.resize((size_t)natom);
  n.resize((size_t)natom);
  int idx = 0;
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at, ++idx) {
    const double* xyz = frm.XYZ(*at);
    double a, b, c;
    ST_SplitXYZ(nax, xyz[0], xyz[1], xyz[2], a, b, c);
    t1[idx] = ST_Wrap(a, Lt1);
    t2[idx] = ST_Wrap(b, Lt2);
    n[idx] = c;
  }
}

/** Number density of mask atoms in the mid-slab slab |n − L_n/2| ≤ w (Å⁻³).
  * Volume is L₁ L₂ (2w). Used for the rhobulk DataSet on the ITIM path;
  * Willard–Chandler ρ_bulk is the laterally averaged coarse-grained field.
  */
static double ST_RhoBulkFromAtoms(std::vector<double> const& ncoord, int natom,
                                  double Lt1, double Lt2, double Ln,
                                  double bulk_halfwidth)
{
  double mid = 0.5 * Ln;
  int nbulk = 0;
  for (int i = 0; i < natom; i++) {
    if (fabs(ncoord[i] - mid) <= bulk_halfwidth)
      nbulk++;
  }
  double vol = Lt1 * Lt2 * 2.0 * bulk_halfwidth;
  if (vol <= 0.0) return 0.0;
  return (double)nbulk / vol;
}

/** Capillary roughness of one instantaneous surface, in Å:
  *   w = √⟨ (h − ⟨h⟩)² ⟩_{i₁,i₂}
  * the RMS of the same field whose DFT gives S(q).
  */
static double ST_RMS(std::vector<double> const& h) {
  if (h.empty()) return ST_NaN();
  double mean = 0.0;
  for (size_t i = 0; i < h.size(); i++)
    mean += h[i];
  mean /= (double)h.size();
  double acc = 0.0;
  for (size_t i = 0; i < h.size(); i++) {
    double d = h[i] - mean;
    acc += d * d;
  }
  return sqrt(acc / (double)h.size());
}

/** Cyclic frequency of DFT bin k on an N-point grid of spacing Δ (Å).
  * Mode index n ∈ (−N/2, N/2]; frequency is n / (N Δ) = n / L_α.
  * Wavevector component is then q_α = 2π n / L_α (Å⁻¹).
  */
static double ST_FftFreq(int k, int n, double d) {
  int p;
  if (k <= (n - 1) / 2)
    p = k;
  else
    p = k - n;
  return (double)p / ((double)n * d);
}

/** |q| = √(q₁² + q₂²) on the N₁×N₂ Fourier mesh, q_α = 2π n_α / L_α. */
static void ST_MakeQGrid(int nx, int ny, double Lx, double Ly, std::vector<double>& q) {
  q.resize((size_t)nx * (size_t)ny);
  double dx = Lx / (double)nx;
  double dy = Ly / (double)ny;
  for (int kx = 0; kx < nx; kx++) {
    double qx = Constants::TWOPI * ST_FftFreq(kx, nx, dx);
    for (int ky = 0; ky < ny; ky++) {
      double qy = Constants::TWOPI * ST_FftFreq(ky, ny, dy);
      q[ST_Idx2(kx, ky, ny)] = sqrt(qx * qx + qy * qy);
    }
  }
}

/// One isotropic |q| shell: all modes with the same rounded q².
struct ST_Shell {
  double q;     ///< Mean |q| of modes in this shell (Å⁻¹)
  double S;     ///< Combined S(q) ≡ ⟨|h_q|²⟩ (Å²)
  double Stop;  ///< Upper-interface S(q)
  double Sbot;  ///< Lower-interface S(q)
  int n;        ///< Number of Fourier modes averaged into this shell
};

/// Sort shells by increasing |q|.
static bool ST_ShellQCmp(ST_Shell const& a, ST_Shell const& b) {
  return a.q < b.q;
}

/** Isotropic shell average of the 2-D spectrum.
  * Modes with the same q² rounded to 10 decimal places (numpy / pandas
  * round(q², 10)) are one shell. The q = 0 (mean-height) mode is dropped.
  * Each shell stores the mean |q| and the mean of S, S_upper, S_lower.
  */
static void ST_ShellAverage(std::vector<double> const& q,
                            std::vector<double> const& combined,
                            std::vector<double> const& top,
                            std::vector<double> const& bot,
                            std::vector<ST_Shell>& shells)
{
  shells.clear();
  struct Acc {
    Acc() : qsum(0.0), S(0.0), Stop(0.0), Sbot(0.0), n(0) {}
    double qsum, S, Stop, Sbot;
    int n;
  };
  std::map<long long, Acc> bins;
  size_t n = q.size();
  for (size_t i = 0; i < n; i++) {
    if (q[i] <= 0.0) continue;
    // Group by rounded q², matching numpy.round(q**2, decimals=10).
    long long key = (long long)floor(q[i] * q[i] * 1.0e10 + 0.5);
    Acc& a = bins[key];
    a.qsum += q[i];
    a.S += combined[i];
    a.Stop += top[i];
    a.Sbot += bot[i];
    a.n++;
  }
  shells.reserve(bins.size());
  for (std::map<long long, Acc>::const_iterator it = bins.begin(); it != bins.end(); ++it) {
    ST_Shell s;
    s.n = it->second.n;
    s.q = it->second.qsum / (double)s.n;
    s.S = it->second.S / (double)s.n;
    s.Stop = it->second.Stop / (double)s.n;
    s.Sbot = it->second.Sbot / (double)s.n;
    shells.push_back(s);
  }
  std::sort(shells.begin(), shells.end(), ST_ShellQCmp);
}

/** Capillary-wave γ from the small-q plateau of q² S(q).
  *   S(q) = k_B T / (A γ q²)                 (κ → 0)
  *   γ    = k_B T / (A ⟨q² S(q)⟩)
  * ⟨…⟩ is an equal-weight mean over shells on [q_min, q_max].
  * A is converted Å² → m²; γ is then ×1000 → mN/m.
  * Needs at least two shells.
  * \return 0 on success, 1 on failure (γ / ⟨q² S⟩ set to NaN).
  */
static int ST_CalcGamma(std::vector<ST_Shell> const& shells, double temperature,
                        double area_A2, double qmin, double qmax,
                        double& gamma, double& plateau)
{
  gamma = ST_NaN();
  plateau = ST_NaN();
  double acc = 0.0;
  int nfit = 0;
  for (size_t i = 0; i < shells.size(); i++) {
    if (shells[i].q >= qmin && shells[i].q <= qmax) {
      acc += shells[i].q * shells[i].q * shells[i].S;
      nfit++;
    }
  }
  if (nfit < 2) return 1;
  plateau = acc / (double)nfit;
  if (plateau <= 0.0) return 1;
  double area_m2 = area_A2 * ST_ANG2_TO_M2;
  gamma = 1000.0 * ST_KB * temperature / (area_m2 * plateau);
  return 0;
}

/** Helfrich inversion of S(q) = k_B T / [A (γ q² + κ q⁴)]:
  *   1/(q² S) = (A / k_B T) (γ + κ q²) = a + b q²
  *   γ         = k_B T a / A            (N/m, then ×1000 → mN/m)
  *   κ / k_B T = b / A                  (A in Å²)
  * q in Å⁻¹, S in Å²; q² S is dimensionless. Equal weight per shell.
  * Needs at least three shells with S > 0.
  * \return 0 on success, 1 on failure (outputs set to NaN).
  */
static int ST_CalcKappa(std::vector<ST_Shell> const& shells, double temperature,
                        double area_A2, double qmin, double qmax,
                        double& gamma, double& kappa_kT)
{
  gamma = ST_NaN();
  kappa_kT = ST_NaN();
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
  int n = 0;
  for (size_t i = 0; i < shells.size(); i++) {
    if (shells[i].q < qmin || shells[i].q > qmax) continue;
    double q2 = shells[i].q * shells[i].q;
    double q2S = q2 * shells[i].S;
    if (q2S <= 0.0 || !ST_Finite(q2S)) continue;
    double x = q2;
    double y = 1.0 / q2S;
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
    n++;
  }
  if (n < 3) return 1;
  double den = (double)n * sxx - sx * sx;
  if (fabs(den) < Constants::SMALL) return 1;
  double a = (sxx * sy - sx * sxy) / den;
  double b = ((double)n * sxy - sx * sy) / den;
  if (a <= 0.0 || !ST_Finite(a) || !ST_Finite(b)) return 1;
  double area_m2 = area_A2 * ST_ANG2_TO_M2;
  gamma = 1000.0 * ST_KB * temperature * a / area_m2;
  kappa_kT = b / area_A2;
  return 0;
}

/** Apparent κ(q)/k_B T at one shell, given a reference γ (mN/m).
  * Invert S = k_B T / [A (γ q² + κ q⁴)] → κ / k_B T = [1/(q² S) − A γ / k_B T] / (A q²).
  * Negative κ is allowed (the linear Helfrich slope may change sign).
  */
static double ST_ShellKappa(double q, double S, double temperature, double area_A2,
                            double gamma_mNm)
{
  if (q <= 0.0 || S <= 0.0 || !ST_Finite(gamma_mNm)) return ST_NaN();
  double q2 = q * q;
  double q2S = q2 * S;
  if (q2S <= 0.0) return ST_NaN();
  double y = 1.0 / q2S;
  double area_m2 = area_A2 * ST_ANG2_TO_M2;
  double a = area_m2 * (gamma_mNm / 1000.0) / (ST_KB * temperature);
  double b = (y - a) / q2;
  return b / area_A2;
}

/** Apparent γ(q) = k_B T / (A q² S(q)) from one shell, in mN/m.
  * Equals the plateau γ only where q² S is flat (κ → 0).
  */
static double ST_ShellGamma(double q, double S, double temperature, double area_A2) {
  double q2S = q * q * S;
  if (q2S <= 0.0) return ST_NaN();
  double area_m2 = area_A2 * ST_ANG2_TO_M2;
  return 1000.0 * ST_KB * temperature / (area_m2 * q2S);
}

/** OLS slope of ln S vs ln q on [q_min, q_max].
  * Pure capillary waves have S(q) ∝ q⁻², so the slope is −2. A slope
  * near −4 is Helfrich-dominated. Sfield selects S / S_upper / S_lower.
  */
static double ST_LogSlope(std::vector<ST_Shell> const& shells, double qmin, double qmax,
                          double ST_Shell::* Sfield)
{
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
  int n = 0;
  for (size_t i = 0; i < shells.size(); i++) {
    double Sval = shells[i].*Sfield;
    if (shells[i].q >= qmin && shells[i].q <= qmax && Sval > 0.0) {
      double x = log(shells[i].q);
      double y = log(Sval);
      sx += x;
      sy += y;
      sxx += x * x;
      sxy += x * y;
      n++;
    }
  }
  if (n < 2) return ST_NaN();
  double den = (double)n * sxx - sx * sx;
  if (fabs(den) < Constants::SMALL) return ST_NaN();
  return ((double)n * sxy - sx * sy) / den;
}

/** Add ds to each non-null DataFile (ASCII / Grace / gnuplot). */
static void ST_AddSetToFiles(DataSet* ds, DataFile* a, DataFile* b, DataFile* c) {
  if (ds == 0) return;
  if (a != 0) a->AddDataSet(ds);
  if (b != 0) b->AddDataSet(ds);
  if (c != 0) c->AddDataSet(ds);
}

/// Write a key/value line; skip non-finite doubles.
static void ST_SummaryD(CpptrajFile* f, const char* key, double v) {
  if (f == 0 || !ST_Finite(v)) return;
  f->Printf("%s %g\n", key, v);
}

static void ST_SummaryI(CpptrajFile* f, const char* key, int v) {
  if (f == 0) return;
  f->Printf("%s %i\n", key, v);
}

static void ST_SummaryS(CpptrajFile* f, const char* key, const char* v) {
  if (f == 0 || v == 0) return;
  f->Printf("%s %s\n", key, v);
}

/** DataFiles are opened only at the end of the run. Fail at Init if the
  * parent directory is missing so a long trajectory is not processed only
  * to lose the write.
  */
static int ST_CheckParentDir(DataFile* df, const char* key) {
  if (df == 0) return 0;
  std::string dir = df->DataFilename().DirPrefix_NoSlash();
  if (dir.empty()) return 0;
  struct stat st;
  if (stat(dir.c_str(), &st) != 0) {
    mprinterr("Error: Directory '%s' for %s does not exist.\n", dir.c_str(), key);
    mprinterr("Error: Create it first, or give a filename in the current directory.\n");
    return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
/// CONSTRUCTOR — Willard–Chandler, nsurf 2, normal z; grid / σ match the Python.
Action_SurfaceTension::Action_SurfaceTension() :
  iface_(WILLARD),
  normal_(AXIS_Z),
  side_(SIDE_UPPER),
  nsurf_(2),
  has_mask2_(false),
  do_upper_(true),
  do_lower_(true),
  qmin_specified_(false),
  temp_(-1.0),
  gridspacing_(2.5),
  dz_(1.0),
  sigma_xy_(2.5),
  sigma_z_(1.5),
  bulk_halfwidth_(5.0),
  threshold_frac_(0.5),
  qmin_(-1.0),
  qmax_(0.174649),
  q_fundamental_(0.0),
  lx_user_(-1.0),
  ly_user_(-1.0),
  lz_user_(-1.0),
  dt_(-1.0),
  blocktime_(-1.0),
  nblock_(0),
  debug_(0),
  n_skip_warn_(0),
  S_(0), S_top_(0), S_bot_(0),
  q2S_(0), q2S_top_(0), q2S_bot_(0),
  gammaq_(0), gammaq_top_(0), gammaq_bot_(0),
  kappaq_(0), kappaq_top_(0), kappaq_bot_(0),
  wtop_(0), wbot_(0), wmean_(0), rhobulk_(0),
  block_gamma_(0), block_kappa_(0), block_wmean_(0), block_wtop_(0), block_wbot_(0),
  summaryFile_(0),
  nx_(0), ny_(0), nz_(0),
  Lt1_ref_(0.0), Lt2_ref_(0.0),
  grid_ready_(false),
  n_frames_(0), n_surfaces_(0), n_skipped_(0), n_blocks_(0),
  block_surface_count_(0), block_frame_count_(0),
  block_w_sum_(0.0), block_wtop_sum_(0.0), block_wbot_sum_(0.0)
{}

// Action_SurfaceTension::Help()
/** ASCII-only (terminals). Formulae live in the class header, not here. */
void Action_SurfaceTension::Help() const {
  mprintf("\t[<name>] <mask> [mask2 <mask2>] temp <T>\n"
          "\t[normal {x|y|z}] [nsurf {1|2}] [side {upper|lower}]\n"
          "\t[interface {willard|itim}]\n"
          "\t[gridspacing <d>] [dz <d> | dnormal <d>]\n"
          "\t[sigmaxy <d>] [sigmaz <d> | sigmanormal <d>]\n"
          "\t[bulkhalfwidth <d>] [threshold <frac>]\n"
          "\t[qmin <q>] [qmax <q>] [lx <Lx>] [ly <Ly>] [lz <Lz>]\n"
          "\t[nblock <frames>] [dt <ps>] [blocktime <ps>]\n"
          "\t[out <file> | spectrumout <file>] [roughout <file>] [blockout <file>]\n"
          "\t[summaryout <file>]\n"
          "\t[spectrumagr <file>] [roughagr <file>] [blockagr <file>]\n"
          "\t[spectrumgnu <file>] [roughgnu <file>] [blockgnu <file>]\n"
          "  Capillary-wave surface tension (mN/m) for a liquid slab.\n"
          "  Lateral box lengths come from the unit cell (NVT). lx/ly/lz only\n"
          "  override noisy box records. Default is two interfaces (vacuum or a\n"
          "  second phase on both sides, e.g. water slab with vacuum at +z and\n"
          "  -z). nsurf 1 is a single interface (side upper|lower). mask2 is a\n"
          "  second atom set for the lower surface (leaflet or second liquid);\n"
          "  the upper surface then comes from <mask> with no mid-box split.\n"
          "  If qmin is omitted it is 2*pi/max(Lt1,Lt2) from the first frame.\n"
          "  blocktime (ps) with dt (analyzed-frame spacing, ps) sets nblock.\n"
          "  out is an alias for spectrumout. Write to the current directory,\n"
          "  e.g. spectrumout spec.dat roughout rough.dat blockout blocks.dat.\n"
          "  A directory in the name must already exist; it is not created.\n");
}

// Action_SurfaceTension::Init()
/** Parse keywords, allocate DataSets, attach optional output files.
  * Spectrum / roughness meshes always exist so writedata can dump them after
  * Print() fills S(q). Files are written only if the matching *out / *agr /
  * *gnu keyword is given. Parent directories are checked here: DataFiles
  * open only after run, so a missing path would otherwise lose the write.
  * out is an alias for spectrumout. agr/gnu force Grace / gnuplot.
  * Contains() is called before getKey* so flags are not already consumed.
  */
Action::RetType Action_SurfaceTension::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  // Optional output files. AddDataFile returns 0 if the keyword is absent.
  // out == spectrumout (usual cpptraj keyword for the main result file).
  bool has_spectrumout = actionArgs.Contains("spectrumout");
  bool has_out = actionArgs.Contains("out");
  std::string specname = actionArgs.GetStringKey("spectrumout");
  std::string outname = actionArgs.GetStringKey("out");
  if (has_spectrumout && has_out && specname != outname) {
    mprinterr("Error: 'out' and 'spectrumout' both specified and differ.\n");
    return Action::ERR;
  }
  if (specname.empty())
    specname = outname;
  bool has_summaryout = actionArgs.Contains("summaryout");
  std::string sumname = actionArgs.GetStringKey("summaryout");
  DataFile* spectrumFile = init.DFL().AddDataFile(specname, actionArgs);
  DataFile* roughFile    = init.DFL().AddDataFile(actionArgs.GetStringKey("roughout"), actionArgs);
  DataFile* blockFile    = init.DFL().AddDataFile(actionArgs.GetStringKey("blockout"), actionArgs);
  DataFile* spectrumAgr  = init.DFL().AddDataFile(actionArgs.GetStringKey("spectrumagr"), actionArgs, DataFile::XMGRACE);
  DataFile* roughAgr     = init.DFL().AddDataFile(actionArgs.GetStringKey("roughagr"), actionArgs, DataFile::XMGRACE);
  DataFile* blockAgr     = init.DFL().AddDataFile(actionArgs.GetStringKey("blockagr"), actionArgs, DataFile::XMGRACE);
  DataFile* spectrumGnu  = init.DFL().AddDataFile(actionArgs.GetStringKey("spectrumgnu"), actionArgs, DataFile::GNUPLOT);
  DataFile* roughGnu     = init.DFL().AddDataFile(actionArgs.GetStringKey("roughgnu"), actionArgs, DataFile::GNUPLOT);
  DataFile* blockGnu     = init.DFL().AddDataFile(actionArgs.GetStringKey("blockgnu"), actionArgs, DataFile::GNUPLOT);
  summaryFile_ = init.DFL().AddCpptrajFile(sumname, "SurfTension summary");
  if (has_summaryout && summaryFile_ == 0) {
    mprinterr("Error: Could not open summaryout '%s'.\n", sumname.c_str());
    return Action::ERR;
  }
  if (ST_CheckParentDir(spectrumFile, "spectrumout") ||
      ST_CheckParentDir(roughFile, "roughout") ||
      ST_CheckParentDir(blockFile, "blockout") ||
      ST_CheckParentDir(spectrumAgr, "spectrumagr") ||
      ST_CheckParentDir(roughAgr, "roughagr") ||
      ST_CheckParentDir(blockAgr, "blockagr") ||
      ST_CheckParentDir(spectrumGnu, "spectrumgnu") ||
      ST_CheckParentDir(roughGnu, "roughgnu") ||
      ST_CheckParentDir(blockGnu, "blockgnu"))
    return Action::ERR;

  temp_ = actionArgs.getKeyDouble("temp", -1.0);
  if (temp_ <= 0.0) {
    mprinterr("Error: Temperature must be specified with 'temp <T>' (T > 0).\n");
    return Action::ERR;
  }
  gridspacing_ = actionArgs.getKeyDouble("gridspacing", 2.5);
  sigma_xy_ = actionArgs.getKeyDouble("sigmaxy", 2.5);
  bulk_halfwidth_ = actionArgs.getKeyDouble("bulkhalfwidth", 5.0);
  threshold_frac_ = actionArgs.getKeyDouble("threshold", 0.5);
  qmin_specified_ = actionArgs.Contains("qmin");
  qmin_ = actionArgs.getKeyDouble("qmin", -1.0);
  qmax_ = actionArgs.getKeyDouble("qmax", 0.174649);
  bool has_lx = actionArgs.Contains("lx");
  bool has_ly = actionArgs.Contains("ly");
  bool has_lz = actionArgs.Contains("lz");
  lx_user_ = actionArgs.getKeyDouble("lx", -1.0);
  ly_user_ = actionArgs.getKeyDouble("ly", -1.0);
  lz_user_ = actionArgs.getKeyDouble("lz", -1.0);
  bool has_nblock = actionArgs.Contains("nblock");
  bool has_blocktime = actionArgs.Contains("blocktime");
  nblock_ = actionArgs.getKeyInt("nblock", 0);
  dt_ = actionArgs.getKeyDouble("dt", -1.0);
  blocktime_ = actionArgs.getKeyDouble("blocktime", -1.0);

  // dz / dnormal and sigmaz / sigmanormal are aliases (normal-axis spacing / sigma).
  bool has_dz = actionArgs.Contains("dz");
  bool has_dnormal = actionArgs.Contains("dnormal");
  double dz_val = actionArgs.getKeyDouble("dz", 1.0);
  double dnormal_val = actionArgs.getKeyDouble("dnormal", 1.0);
  if (has_dz && has_dnormal && fabs(dz_val - dnormal_val) > 1.0e-12) {
    mprinterr("Error: dz and dnormal both specified and differ.\n");
    return Action::ERR;
  }
  dz_ = has_dnormal ? dnormal_val : dz_val;

  bool has_sigmaz = actionArgs.Contains("sigmaz");
  bool has_sigman = actionArgs.Contains("sigmanormal");
  double sigmaz_val = actionArgs.getKeyDouble("sigmaz", 1.5);
  double sigman_val = actionArgs.getKeyDouble("sigmanormal", 1.5);
  if (has_sigmaz && has_sigman && fabs(sigmaz_val - sigman_val) > 1.0e-12) {
    mprinterr("Error: sigmaz and sigmanormal both specified and differ.\n");
    return Action::ERR;
  }
  sigma_z_ = has_sigman ? sigman_val : sigmaz_val;

  std::string nstr = actionArgs.GetStringKey("normal");
  if (nstr.empty() && actionArgs.Contains("normal")) {
    mprinterr("Error: 'normal' requires x, y, or z.\n");
    return Action::ERR;
  }
  if (nstr.empty())
    nstr = "z";
  if (nstr == "x" || nstr == "X")
    normal_ = AXIS_X;
  else if (nstr == "y" || nstr == "Y")
    normal_ = AXIS_Y;
  else if (nstr == "z" || nstr == "Z")
    normal_ = AXIS_Z;
  else {
    mprinterr("Error: normal must be 'x', 'y', or 'z'.\n");
    return Action::ERR;
  }

  std::string ifacestr = actionArgs.GetStringKey("interface");
  if (ifacestr.empty() && actionArgs.Contains("interface")) {
    mprinterr("Error: interface must be 'willard' or 'itim'.\n");
    return Action::ERR;
  }
  if (ifacestr.empty())
    ifacestr = "willard";
  if (ifacestr == "willard" || ifacestr == "wc")
    iface_ = WILLARD;
  else if (ifacestr == "itim" || ifacestr == "minmax")
    iface_ = ITIM;
  else {
    mprinterr("Error: interface must be 'willard' or 'itim'.\n");
    return Action::ERR;
  }

  nsurf_ = actionArgs.getKeyInt("nsurf", 2);
  if (nsurf_ != 1 && nsurf_ != 2) {
    mprinterr("Error: nsurf must be 1 or 2.\n");
    return Action::ERR;
  }
  bool has_side = actionArgs.Contains("side");
  std::string sidestr = actionArgs.GetStringKey("side");
  if (sidestr.empty() && has_side) {
    mprinterr("Error: 'side' requires upper or lower.\n");
    return Action::ERR;
  }
  if (sidestr.empty())
    sidestr = "upper";
  if (sidestr == "upper" || sidestr == "top")
    side_ = SIDE_UPPER;
  else if (sidestr == "lower" || sidestr == "bot" || sidestr == "bottom")
    side_ = SIDE_LOWER;
  else {
    mprinterr("Error: side must be 'upper' or 'lower'.\n");
    return Action::ERR;
  }
  if (nsurf_ == 2) {
    do_upper_ = true;
    do_lower_ = true;
    if (has_side)
      mprintf("Warning: 'side' is ignored when nsurf is 2.\n");
  } else {
    do_upper_ = (side_ == SIDE_UPPER);
    do_lower_ = (side_ == SIDE_LOWER);
  }

  std::string mask2exp = actionArgs.GetStringKey("mask2");
  has_mask2_ = !mask2exp.empty();
  if (has_mask2_) {
    if (nsurf_ != 2) {
      mprinterr("Error: mask2 requires nsurf 2 (upper from <mask>, lower from mask2).\n");
      return Action::ERR;
    }
    if (Mask2_.SetMaskString(mask2exp)) return Action::ERR;
  }

  if (has_lx && lx_user_ <= 0.0) {
    mprinterr("Error: lx must be > 0.\n");
    return Action::ERR;
  }
  if (has_ly && ly_user_ <= 0.0) {
    mprinterr("Error: ly must be > 0.\n");
    return Action::ERR;
  }
  if (has_lz && lz_user_ <= 0.0) {
    mprinterr("Error: lz must be > 0.\n");
    return Action::ERR;
  }

  if (gridspacing_ <= 0.0) {
    mprinterr("Error: gridspacing must be > 0.\n");
    return Action::ERR;
  }
  if (iface_ == WILLARD) {
    if (dz_ <= 0.0) {
      mprinterr("Error: dz (or dnormal) must be > 0.\n");
      return Action::ERR;
    }
    if (sigma_xy_ < 0.0 || sigma_z_ < 0.0) {
      mprinterr("Error: sigmaxy and sigmaz (or sigmanormal) must be >= 0.\n");
      return Action::ERR;
    }
    if (threshold_frac_ <= 0.0 || threshold_frac_ >= 1.0) {
      mprinterr("Error: threshold must be between 0 and 1.\n");
      return Action::ERR;
    }
  }
  if (bulk_halfwidth_ <= 0.0) {
    mprinterr("Error: bulkhalfwidth must be > 0.\n");
    return Action::ERR;
  }
  if (qmax_ <= 0.0) {
    mprinterr("Error: qmax must be > 0.\n");
    return Action::ERR;
  }
  if (qmin_specified_) {
    if (qmin_ <= 0.0) {
      mprinterr("Error: qmin must be > 0.\n");
      return Action::ERR;
    }
    if (qmax_ <= qmin_) {
      mprinterr("Error: qmax must be greater than qmin.\n");
      return Action::ERR;
    }
  }
  if (has_blocktime) {
    if (dt_ <= 0.0) {
      mprinterr("Error: 'blocktime' requires 'dt <ps>' (time between analyzed frames).\n");
      return Action::ERR;
    }
    if (blocktime_ <= 0.0) {
      mprinterr("Error: blocktime must be > 0.\n");
      return Action::ERR;
    }
    int nfromtime = ST_IRound(blocktime_ / dt_);
    if (nfromtime < 1) {
      mprinterr("Error: blocktime / dt is less than one frame.\n");
      return Action::ERR;
    }
    if (has_nblock && nblock_ != nfromtime) {
      mprinterr("Error: nblock (%i) does not match blocktime/dt (%i frames).\n",
                nblock_, nfromtime);
      return Action::ERR;
    }
    nblock_ = nfromtime;
  }
  if (nblock_ < 0) {
    mprinterr("Error: nblock must be >= 0.\n");
    return Action::ERR;
  }
  if ((blockFile != 0 || blockAgr != 0 || blockGnu != 0) && nblock_ < 1) {
    mprinterr("Error: 'blockout'/'blockagr'/'blockgnu' require 'nblock' or 'blocktime'.\n");
    return Action::ERR;
  }

  std::string maskexp = actionArgs.GetMaskNext();
  if (maskexp.empty()) {
    mprinterr("Error: No atom mask given.\n");
    return Action::ERR;
  }
  if (Mask_.SetMaskString(maskexp)) return Action::ERR;

  // Dataset name is the leftover argument (cpptraj convention). Reuse the
  // generated name so all aspects share one set family (ST_00000, …).
  std::string dsname = actionArgs.GetStringNext();
  S_ = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "S", MetaData::NOT_TS), "ST");
  if (S_ == 0) return Action::ERR;
  dsname = S_->Meta().Name();
  S_top_     = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "Stop", MetaData::NOT_TS), "ST");
  S_bot_     = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "Sbot", MetaData::NOT_TS), "ST");
  q2S_       = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "q2S", MetaData::NOT_TS), "ST");
  q2S_top_   = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "q2Stop", MetaData::NOT_TS), "ST");
  q2S_bot_   = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "q2Sbot", MetaData::NOT_TS), "ST");
  gammaq_    = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "gammaq", MetaData::NOT_TS), "ST");
  gammaq_top_= (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "gammaqtop", MetaData::NOT_TS), "ST");
  gammaq_bot_= (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "gammaqbot", MetaData::NOT_TS), "ST");
  kappaq_    = (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "kappaq", MetaData::NOT_TS), "ST");
  kappaq_top_= (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "kappaqtop", MetaData::NOT_TS), "ST");
  kappaq_bot_= (DataSet_Mesh*)init.DSL().AddSet(DataSet::XYMESH, MetaData(dsname, "kappaqbot", MetaData::NOT_TS), "ST");
  if (S_==0 || S_top_==0 || S_bot_==0 || q2S_==0 || q2S_top_==0 || q2S_bot_==0 ||
      gammaq_==0 || gammaq_top_==0 || gammaq_bot_==0 ||
      kappaq_==0 || kappaq_top_==0 || kappaq_bot_==0)
    return Action::ERR;
  // Grace/gnuplot xlabel comes from Dim(0). Mesh X values are q, not a uniform grid.
  Dimension qdim(0.0, 1.0, "q (Ang^-1)");
  S_->SetDim(Dimension::X, qdim);
  S_top_->SetDim(Dimension::X, qdim);
  S_bot_->SetDim(Dimension::X, qdim);
  q2S_->SetDim(Dimension::X, qdim);
  q2S_top_->SetDim(Dimension::X, qdim);
  q2S_bot_->SetDim(Dimension::X, qdim);
  gammaq_->SetDim(Dimension::X, qdim);
  gammaq_top_->SetDim(Dimension::X, qdim);
  gammaq_bot_->SetDim(Dimension::X, qdim);
  kappaq_->SetDim(Dimension::X, qdim);
  kappaq_top_->SetDim(Dimension::X, qdim);
  kappaq_bot_->SetDim(Dimension::X, qdim);
  // Print() fills these from the accumulated spectra; they are not time series.
  // MPI: SyncAction reduces |h_q|²; skip DataSet concat of empty meshes.
# ifdef MPI
  S_->SetNeedsSync(false); S_top_->SetNeedsSync(false); S_bot_->SetNeedsSync(false);
  q2S_->SetNeedsSync(false); q2S_top_->SetNeedsSync(false); q2S_bot_->SetNeedsSync(false);
  gammaq_->SetNeedsSync(false); gammaq_top_->SetNeedsSync(false); gammaq_bot_->SetNeedsSync(false);
  kappaq_->SetNeedsSync(false); kappaq_top_->SetNeedsSync(false); kappaq_bot_->SetNeedsSync(false);
# endif
  ST_AddSetToFiles(S_,          spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(S_top_,      spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(S_bot_,      spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(q2S_,        spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(q2S_top_,    spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(q2S_bot_,    spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(gammaq_,     spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(gammaq_top_, spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(gammaq_bot_, spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(kappaq_,     spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(kappaq_top_, spectrumFile, spectrumAgr, spectrumGnu);
  ST_AddSetToFiles(kappaq_bot_, spectrumFile, spectrumAgr, spectrumGnu);

  wtop_   = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "wtop"), "ST");
  wbot_   = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "wbot"), "ST");
  wmean_  = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "wmean"), "ST");
  rhobulk_= init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "rhobulk"), "ST");
  if (wtop_==0 || wbot_==0 || wmean_==0 || rhobulk_==0) return Action::ERR;
  ST_AddSetToFiles(wtop_,    roughFile, roughAgr, roughGnu);
  ST_AddSetToFiles(wbot_,    roughFile, roughAgr, roughGnu);
  ST_AddSetToFiles(wmean_,   roughFile, roughAgr, roughGnu);
  ST_AddSetToFiles(rhobulk_, roughFile, roughAgr, roughGnu);

  if (nblock_ > 0) {
    // Block DataSets are indexed by completed block, not by frame.
    block_gamma_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "bgamma"), "ST");
    block_kappa_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "bkappa"), "ST");
    block_wmean_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "bwmean"), "ST");
    block_wtop_  = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "bwtop"), "ST");
    block_wbot_  = init.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "bwbot"), "ST");
    if (block_gamma_==0 || block_kappa_==0 || block_wmean_==0 || block_wtop_==0 || block_wbot_==0)
      return Action::ERR;
    ST_AddSetToFiles(block_gamma_, blockFile, blockAgr, blockGnu);
    ST_AddSetToFiles(block_kappa_, blockFile, blockAgr, blockGnu);
    ST_AddSetToFiles(block_wmean_, blockFile, blockAgr, blockGnu);
    ST_AddSetToFiles(block_wtop_,  blockFile, blockAgr, blockGnu);
    ST_AddSetToFiles(block_wbot_,  blockFile, blockAgr, blockGnu);
  }

  mprintf("    SURFTENSION: Capillary-wave surface tension.\n");
  mprintf("\tMask: '%s'\n", Mask_.MaskString());
  if (has_mask2_)
    mprintf("\tMask2 (lower surface): '%s'\n", Mask2_.MaskString());
  {
    const char* t1 = 0;
    const char* t2 = 0;
    ST_PlaneAxes((int)normal_, t1, t2);
    mprintf("\tSlab normal: %s (interface plane %s-%s)\n",
            ST_AxisName((int)normal_), t1, t2);
  }
  if (nsurf_ == 2)
    mprintf("\tSurfaces: 2 (upper and lower; vacuum/second phase on both sides).\n");
  else
    mprintf("\tSurfaces: 1 (%s only).\n", (side_ == SIDE_UPPER) ? "upper" : "lower");
  if (iface_ == WILLARD)
    mprintf("\tInterface: Willard-Chandler density isosurface.\n");
  else if (has_mask2_)
    mprintf("\tInterface: ITIM (upper = max of mask, lower = min of mask2).\n");
  else
    mprintf("\tInterface: ITIM min/max (per-column, split at mid-box along %s).\n",
            ST_AxisName((int)normal_));
  mprintf("\tTemperature= %g K\n", temp_);
  mprintf("\tLateral box lengths: from the unit cell");
  if (lx_user_ > 0.0 || ly_user_ > 0.0 || lz_user_ > 0.0)
    mprintf(" (some Cartesian lengths overridden)");
  mprintf(".\n");
  mprintf("\tInterface-plane grid spacing= %g Ang\n", gridspacing_);
  if (iface_ == WILLARD) {
    mprintf("\tNormal-axis bin spacing (dz)= %g Ang\n", dz_);
    mprintf("\tGaussian sigma in plane= %g Ang, along normal= %g Ang\n",
            sigma_xy_, sigma_z_);
    mprintf("\tBulk half-width= %g Ang, threshold fraction= %g\n",
            bulk_halfwidth_, threshold_frac_);
  } else {
    mprintf("\tBulk half-width (rho_bulk only)= %g Ang\n", bulk_halfwidth_);
  }
  if (qmin_specified_)
    mprintf("\tFit q range= %g to %g Ang^-1\n", qmin_, qmax_);
  else
    mprintf("\tFit qmin: 2*pi/max(Lt1,Lt2) from the first frame; qmax= %g Ang^-1\n",
            qmax_);
  mprintf("\tHelfrich kappa: linear fit of 1/(q^2 S) vs q^2 on that window.\n");
  mprintf("\tHeight-field 2-D FFT: PubFFT (numpy fft2 / (nx*ny)).\n");
  if (lx_user_ > 0.0)
    mprintf("\tUsing fixed Lx= %g Ang (%s)\n", lx_user_,
            (normal_ == AXIS_X) ? "normal" : "lateral");
  if (ly_user_ > 0.0)
    mprintf("\tUsing fixed Ly= %g Ang (%s)\n", ly_user_,
            (normal_ == AXIS_Y) ? "normal" : "lateral");
  if (lz_user_ > 0.0)
    mprintf("\tUsing fixed Lz= %g Ang (%s)\n", lz_user_,
            (normal_ == AXIS_Z) ? "normal" : "lateral");
  if (nblock_ > 0) {
    mprintf("\tBlock averaging every %i analyzed frames", nblock_);
    if (blocktime_ > 0.0 && dt_ > 0.0)
      mprintf(" (blocktime %g ps, dt %g ps)", blocktime_, dt_);
#   ifdef MPI
    if (trajComm_.Size() > 1)
      mprintf(" (per rank)");
#   endif
    mprintf(".\n");
  } else
    mprintf("\tBlock averaging disabled.\n");
  mprintf("\tSpectrum DataSets: %s %s %s %s %s %s %s %s %s %s %s %s\n",
          S_->legend(), S_top_->legend(), S_bot_->legend(),
          q2S_->legend(), q2S_top_->legend(), q2S_bot_->legend(),
          gammaq_->legend(), gammaq_top_->legend(), gammaq_bot_->legend(),
          kappaq_->legend(), kappaq_top_->legend(), kappaq_bot_->legend());
  if (spectrumFile != 0)
    mprintf("\tSpectrum output to '%s' (%s)\n", spectrumFile->DataFilename().full(),
            spectrumFile->FormatString());
  if (spectrumAgr != 0)
    mprintf("\tSpectrum Grace output to '%s'\n", spectrumAgr->DataFilename().full());
  if (spectrumGnu != 0)
    mprintf("\tSpectrum gnuplot output to '%s'\n", spectrumGnu->DataFilename().full());
  if (roughFile != 0)
    mprintf("\tRoughness output to '%s' (%s)\n", roughFile->DataFilename().full(),
            roughFile->FormatString());
  if (roughAgr != 0)
    mprintf("\tRoughness Grace output to '%s'\n", roughAgr->DataFilename().full());
  if (roughGnu != 0)
    mprintf("\tRoughness gnuplot output to '%s'\n", roughGnu->DataFilename().full());
  if (blockFile != 0)
    mprintf("\tBlock output to '%s' (%s)\n", blockFile->DataFilename().full(),
            blockFile->FormatString());
  if (blockAgr != 0)
    mprintf("\tBlock Grace output to '%s'\n", blockAgr->DataFilename().full());
  if (blockGnu != 0)
    mprintf("\tBlock gnuplot output to '%s'\n", blockGnu->DataFilename().full());
  if (summaryFile_ != 0)
    mprintf("\tSummary output to '%s'\n", summaryFile_->Filename().full());
# ifdef _OPENMP
  {
    int nthreads = 1;
#   pragma omp parallel
    {
#     pragma omp master
      nthreads = omp_get_num_threads();
    }
    if (nthreads > 1) {
      if (iface_ == WILLARD)
        mprintf("\tOpenMP: 3-D Gaussian filter parallelized with %i threads.\n", nthreads);
    }
  }
# endif
# ifdef MPI
  if (trajComm_.Size() > 1)
    mprintf("\tMPI: |h_q|^2 spectra reduced to master with one packed SUM.\n");
# endif
  return Action::OK;
}

// Action_SurfaceTension::Setup()
/** Require an orthorhombic unit cell (L₁, L₂, L_n) and bind the mask(s).
  * Non-X-aligned boxes are warned: laterals are wrapped independently.
  */
Action::RetType Action_SurfaceTension::Setup(ActionSetup& setup)
{
  if (!setup.CoordInfo().HasBox()) {
    mprintf("Warning: No unit cell; surface tension cannot be calculated for '%s'\n",
            setup.Top().c_str());
    return Action::SKIP;
  }
  if (!setup.CoordInfo().TrajBox().Is_X_Aligned_Ortho()) {
    mprintf("Warning: Box is not X-aligned orthorhombic; wrapping Cartesian\n"
            "Warning:   coordinates independently may not be correct.\n");
  }
  if (setup.Top().SetupIntegerMask(Mask_)) return Action::ERR;
  Mask_.MaskInfo();
  if (Mask_.None()) {
    mprintf("Warning: Mask '%s' selects no atoms.\n", Mask_.MaskString());
    return Action::SKIP;
  }
  if (has_mask2_) {
    if (setup.Top().SetupIntegerMask(Mask2_)) return Action::ERR;
    Mask2_.MaskInfo();
    if (Mask2_.None()) {
      mprintf("Warning: Mask2 '%s' selects no atoms.\n", Mask2_.MaskString());
      return Action::SKIP;
    }
  }
  return Action::OK;
}

// Action_SurfaceTension::AllocateGrid()
/** First good frame: freeze N₁, N₂, L₁, L₂ and (Willard) N_n.
  * Allocates ρ, h, |h_q|² accumulators, PubFFT plans, and q_α = 2π n_α / L_α.
  * If q_min was omitted it is set to 2π / max(L₁, L₂).
  * \return 0 OK, 1 fatal (FFT setup or q_max ≤ q_min).
  */
int Action_SurfaceTension::AllocateGrid(int nx, int ny, int nz) {
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  size_t n3 = (size_t)nx * (size_t)ny * (size_t)nz;
  size_t n2 = (size_t)nx * (size_t)ny;
  if (nz > 0)
    density_.assign(n3, 0.0);
  else
    density_.clear();
  if (nz > 0 && has_mask2_)
    density2_.assign(n3, 0.0);
  else
    density2_.clear();
  h_upper_.assign(n2, 0.0);
  h_lower_.assign(n2, 0.0);
  total_power_.assign(n2, 0.0);
  top_power_.assign(n2, 0.0);
  bottom_power_.assign(n2, 0.0);
  block_power_.assign(n2, 0.0);
  if (fft_n1_.SetupFFTforN(nx_) || fft_n2_.SetupFFTforN(ny_)) {
    mprinterr("Error: Could not set up 2-D FFT (nx = %i, ny = %i).\n", nx_, ny_);
    return 1;
  }
  fft_grid_.Allocate(nx_ * ny_);
  fft_row_.Allocate(ny_);
  fft_col_.Allocate(nx_);
  ST_MakeQGrid(nx_, ny_, Lt1_ref_, Lt2_ref_, q_grid_);
  double qmin_acc = 0.0;
  bool haveq = false;
  for (size_t i = 0; i < q_grid_.size(); i++) {
    if (q_grid_[i] > 0.0 && (!haveq || q_grid_[i] < qmin_acc)) {
      qmin_acc = q_grid_[i];
      haveq = true;
    }
  }
  const char* t1 = 0;
  const char* t2 = 0;
  ST_PlaneAxes((int)normal_, t1, t2);
  mprintf("\tInterface grid = %i x %i", nx_, ny_);
  if (nz_ > 0)
    mprintf(" x %i", nz_);
  mprintf(" bins along %s, %s", t1, t2);
  if (nz_ > 0)
    mprintf(", %s", ST_AxisName((int)normal_));
  mprintf("\n");
  mprintf("\tL(%s) = %g Ang, L(%s) = %g Ang\n", t1, Lt1_ref_, t2, Lt2_ref_);
  {
    double q1 = Constants::TWOPI / Lt1_ref_;
    double q2 = Constants::TWOPI / Lt2_ref_;
    q_fundamental_ = (q1 < q2) ? q1 : q2;
    mprintf("\tFundamental q: 2*pi/L(%s) = %g, 2*pi/L(%s) = %g Ang^-1\n",
            t1, q1, t2, q2);
    mprintf("\tSmallest fundamental |q| = %g Ang^-1\n", q_fundamental_);
    if (!qmin_specified_) {
      qmin_ = q_fundamental_;
      mprintf("\tFit qmin set from the box = %g Ang^-1\n", qmin_);
    }
    if (qmax_ <= qmin_) {
      mprinterr("Error: qmax (%g) must be greater than qmin (%g).\n", qmax_, qmin_);
      return 1;
    }
  }
  if (haveq)
    mprintf("\tSmallest accessible q (FFT grid) = %g Ang^-1\n", qmin_acc);
  grid_ready_ = true;
  return 0;
}

/** Cap skip-frame warnings at 5, then one “further messages suppressed.” */
void Action_SurfaceTension::SkipWarn(const char* msg)
{
  const int maxw = 5;
  if (n_skip_warn_ < maxw)
    mprintf("Warning: %s\n", msg);
  else if (n_skip_warn_ == maxw)
    mprintf("Warning: Further skip-frame messages suppressed.\n");
  n_skip_warn_++;
}

// Action_SurfaceTension::WillardHeights()
/** Willard–Chandler heights from one atom set into the given density buffer.
  * Histogram atoms onto N₁×N₂×N_n, convert to number density (Å⁻³), apply
  * the periodic Gaussian, take ρ_bulk as the mid-slab lateral mean, then
  * find the ±n crossings of θ ρ_bulk. Writes h_upper_ / h_lower_.
  * \return 0 OK, 1 skip (no bulk bins, ρ_bulk ≤ 0, or a column has no crossing).
  */
int Action_SurfaceTension::WillardHeights(std::vector<double> const& t1,
                                          std::vector<double> const& t2,
                                          std::vector<double> const& ncoord,
                                          int natom,
                                          double Lt1, double Lt2, double Ln,
                                          std::vector<double>& density,
                                          bool want_u, bool want_l,
                                          double& rho_bulk)
{
  double dx = Lt1 / (double)nx_;
  double dy = Lt2 / (double)ny_;
  double dz = Ln / (double)nz_;
  double voxel = dx * dy * dz;
  std::fill(density.begin(), density.end(), 0.0);
  for (int i = 0; i < natom; i++) {
    int ix = (int)floor(t1[i] / dx);
    int iy = (int)floor(t2[i] / dy);
    int iz = (int)floor(ncoord[i] / dz);
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (iz < 0) iz = 0;
    if (ix >= nx_) ix = nx_ - 1;
    if (iy >= ny_) iy = ny_ - 1;
    if (iz >= nz_) iz = nz_ - 1;
    density[ST_Idx3(ix, iy, iz, ny_, nz_)] += 1.0;
  }
  if (voxel > 0.0) {
    for (size_t i = 0; i < density.size(); i++)
      density[i] /= voxel;
  }
  ST_GaussianFilter3D(density, nx_, ny_, nz_,
                      sigma_xy_ / dx, sigma_xy_ / dy, sigma_z_ / dz);

  std::vector<double> n_grid((size_t)nz_);
  for (int iz = 0; iz < nz_; iz++)
    n_grid[iz] = ((double)iz + 0.5) * dz;
  double slab_center = 0.5 * Ln;
  int center_index = 0;
  double best = fabs(n_grid[0] - slab_center);
  for (int iz = 1; iz < nz_; iz++) {
    double d = fabs(n_grid[iz] - slab_center);
    if (d < best) {
      best = d;
      center_index = iz;
    }
  }

  std::vector<double> rho_n((size_t)nz_, 0.0);
  double nxy = (double)(nx_ * ny_);
  for (int ix = 0; ix < nx_; ix++) {
    for (int iy = 0; iy < ny_; iy++) {
      for (int iz = 0; iz < nz_; iz++)
        rho_n[iz] += density[ST_Idx3(ix, iy, iz, ny_, nz_)];
    }
  }
  for (int iz = 0; iz < nz_; iz++)
    rho_n[iz] /= nxy;

  double rho_bulk_acc = 0.0;
  int nbulk = 0;
  for (int iz = 0; iz < nz_; iz++) {
    if (fabs(n_grid[iz] - slab_center) <= bulk_halfwidth_) {
      rho_bulk_acc += rho_n[iz];
      nbulk++;
    }
  }
  if (nbulk < 1) {
    SkipWarn("No bins along the normal fall within the bulk region; skipping frame.");
    return 1;
  }
  rho_bulk = rho_bulk_acc / (double)nbulk;
  if (rho_bulk <= 0.0) {
    SkipWarn("Bulk density is non-positive; skipping frame.");
    return 1;
  }
  double threshold = threshold_frac_ * rho_bulk;

  std::vector<double> col((size_t)nz_);
  for (int ix = 0; ix < nx_; ix++) {
    for (int iy = 0; iy < ny_; iy++) {
      for (int iz = 0; iz < nz_; iz++)
        col[iz] = density[ST_Idx3(ix, iy, iz, ny_, nz_)];
      if (want_u) {
        double hu = ST_FindCrossing(n_grid, col, threshold, center_index, true);
        if (!ST_Finite(hu)) {
          SkipWarn("Could not identify a local interface; skipping frame.");
          return 1;
        }
        h_upper_[ST_Idx2(ix, iy, ny_)] = hu;
      }
      if (want_l) {
        double hl = ST_FindCrossing(n_grid, col, threshold, center_index, false);
        if (!ST_Finite(hl)) {
          SkipWarn("Could not identify a local interface; skipping frame.");
          return 1;
        }
        h_lower_[ST_Idx2(ix, iy, ny_)] = hl;
      }
    }
  }
  return 0;
}

// Action_SurfaceTension::HeightPower()
/** Discrete Fourier transform of h − ⟨h⟩ (row–column 1-D PubFFTs).
  *   h_q = (1/(N₁ N₂)) Σ (h − ⟨h⟩) exp(−i q · r)
  * implemented as unnormalized FFTW / numpy fft2, then ÷ (N₁ N₂).
  * Power is |h_q|² (Å²). fft_n2_ along n₂; fft_n1_ along n₁.
  */
void Action_SurfaceTension::HeightPower(std::vector<double> const& h,
                                        std::vector<double>& power)
{
  int nx = nx_;
  int ny = ny_;
  int nxy = nx * ny;
  power.assign((size_t)nxy, 0.0);
  if (nxy == 0) return;

  double mean = 0.0;
  for (int i = 0; i < nxy; i++)
    mean += h[i];
  mean /= (double)nxy;

  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      int idx = ST_Idx2(ix, iy, ny);
      fft_grid_[2 * idx]     = h[idx] - mean;
      fft_grid_[2 * idx + 1] = 0.0;
    }
  }

  for (int ix = 0; ix < nx; ix++) {
    int base = ix * ny;
    for (int iy = 0; iy < ny; iy++) {
      fft_row_[2 * iy]     = fft_grid_[2 * (base + iy)];
      fft_row_[2 * iy + 1] = fft_grid_[2 * (base + iy) + 1];
    }
    fft_n2_.Forward(fft_row_);
    for (int iy = 0; iy < ny; iy++) {
      fft_grid_[2 * (base + iy)]     = fft_row_[2 * iy];
      fft_grid_[2 * (base + iy) + 1] = fft_row_[2 * iy + 1];
    }
  }

  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      int idx = ST_Idx2(ix, iy, ny);
      fft_col_[2 * ix]     = fft_grid_[2 * idx];
      fft_col_[2 * ix + 1] = fft_grid_[2 * idx + 1];
    }
    fft_n1_.Forward(fft_col_);
    for (int ix = 0; ix < nx; ix++) {
      int idx = ST_Idx2(ix, iy, ny);
      fft_grid_[2 * idx]     = fft_col_[2 * ix];
      fft_grid_[2 * idx + 1] = fft_col_[2 * ix + 1];
    }
  }

  double Ninv = 1.0 / (double)nxy;
  for (int k = 0; k < nxy; k++) {
    double re = fft_grid_[2 * k] * Ninv;
    double im = fft_grid_[2 * k + 1] * Ninv;
    power[k] = re * re + im * im;
  }
}

// Action_SurfaceTension::ProcessFrame()
/** One NVT frame: permute → wrap / recenter → h(t1,t2) → w and |h_q|².
  * First good frame freezes N₁, N₂, L₁, L₂ (and N_n for Willard). Later
  * frames must keep the same lateral box (NVT). Combined power is summed
  * over the surfaces actually used (nsurf 1 or 2).
  * \return 0 OK, 1 skip frame, 2 fatal (grid or L₁, L₂ changed).
  */
int Action_SurfaceTension::ProcessFrame(Frame const& frm, double Lx, double Ly, double Lz)
{
  int nax = (int)normal_;
  double Lt1, Lt2, Ln;
  ST_SplitBox(nax, Lx, Ly, Lz, Lt1, Lt2, Ln);

  std::vector<double> t1, t2, n;
  ST_SplitAtomCoords(frm, Mask_, nax, Lx, Ly, Lz, t1, t2, n);
  int natom = Mask_.Nselected();

  std::vector<double> t1b, t2b, nb;
  int natom2 = 0;
  if (has_mask2_) {
    ST_SplitAtomCoords(frm, Mask2_, nax, Lx, Ly, Lz, t1b, t2b, nb);
    natom2 = Mask2_.Nselected();
    double slab_center = ST_CircularSlabCenter(n, natom, &nb, natom2, Ln);
    ST_ApplyCircularRecenter(n, Ln, slab_center);
    ST_ApplyCircularRecenter(nb, Ln, slab_center);
  } else {
    ST_CircularRecenter(n, Ln);
  }

  int nx = std::max(8, ST_IRound(Lt1 / gridspacing_));
  int ny = std::max(8, ST_IRound(Lt2 / gridspacing_));
  int nz = 0;
  if (iface_ == WILLARD)
    nz = std::max(32, ST_IRound(Ln / dz_));

  if (!grid_ready_) {
    Lt1_ref_ = Lt1;
    Lt2_ref_ = Lt2;
    if (AllocateGrid(nx, ny, nz)) return 2;
  } else {
    if (nx != nx_ || ny != ny_) {
      mprinterr("Error: Interface grid dimensions changed during the trajectory.\n");
      return 2;
    }
    if (fabs(Lt1 - Lt1_ref_) > 1.0e-5 || fabs(Lt2 - Lt2_ref_) > 1.0e-5) {
      mprinterr("Error: Lateral box dimensions changed. surftension assumes a fixed interface plane (NVT).\n");
      return 2;
    }
    if (iface_ == WILLARD && nz != nz_)
      nz = nz_;
  }

  double rho_bulk = 0.0;
  if (iface_ == ITIM) {
    if (has_mask2_) {
      rho_bulk = ST_RhoBulkFromAtoms(n, natom, Lt1, Lt2, Ln, bulk_halfwidth_);
      if (!ST_ItimTwoMasks(t1, t2, n, natom, t1b, t2b, nb, natom2,
                           Lt1, Lt2, nx_, ny_, h_upper_, h_lower_)) {
        SkipWarn("Empty ITIM column; skipping frame.");
        return 1;
      }
    } else {
      rho_bulk = ST_RhoBulkFromAtoms(n, natom, Lt1, Lt2, Ln, bulk_halfwidth_);
      if (!ST_ItimMinMax(t1, t2, n, natom, Lt1, Lt2, Ln, nx_, ny_,
                         do_upper_, do_lower_, h_upper_, h_lower_)) {
        SkipWarn("Empty ITIM column; skipping frame.");
        return 1;
      }
    }
  } else if (has_mask2_) {
    double rho1 = 0.0, rho2 = 0.0;
    if (WillardHeights(t1, t2, n, natom, Lt1, Lt2, Ln, density_, true, false, rho1))
      return 1;
    if (WillardHeights(t1b, t2b, nb, natom2, Lt1, Lt2, Ln, density2_, false, true, rho2))
      return 1;
    rho_bulk = 0.5 * (rho1 + rho2);
  } else {
    if (WillardHeights(t1, t2, n, natom, Lt1, Lt2, Ln, density_,
                       do_upper_, do_lower_, rho_bulk))
      return 1;
  }

  double w_top = do_upper_ ? ST_RMS(h_upper_) : ST_NaN();
  double w_bot = do_lower_ ? ST_RMS(h_lower_) : ST_NaN();
  double w_mean;
  if (do_upper_ && do_lower_)
    w_mean = 0.5 * (w_top + w_bot);
  else
    w_mean = do_upper_ ? w_top : w_bot;

  std::vector<double> p_upper, p_lower;
  int nsurf_this = 0;
  if (do_upper_) {
    HeightPower(h_upper_, p_upper);
    nsurf_this++;
  } else {
    p_upper.assign(top_power_.size(), 0.0);
  }
  if (do_lower_) {
    HeightPower(h_lower_, p_lower);
    nsurf_this++;
  } else {
    p_lower.assign(bottom_power_.size(), 0.0);
  }

  for (size_t i = 0; i < top_power_.size(); i++) {
    if (do_upper_) top_power_[i] += p_upper[i];
    if (do_lower_) bottom_power_[i] += p_lower[i];
    total_power_[i] += p_upper[i] + p_lower[i];
    block_power_[i] += p_upper[i] + p_lower[i];
  }

  n_frames_++;
  n_surfaces_ += nsurf_this;
  block_surface_count_ += nsurf_this;
  block_frame_count_++;
  block_w_sum_ += w_mean;
  block_wtop_sum_ += w_top;
  block_wbot_sum_ += w_bot;

  int fidx = n_frames_ - 1;
  wtop_->Add(fidx, &w_top);
  wbot_->Add(fidx, &w_bot);
  wmean_->Add(fidx, &w_mean);
  rhobulk_->Add(fidx, &rho_bulk);

  if (nblock_ > 0 && (n_frames_ % nblock_) == 0) {
    if (FinishBlock()) return 2;
  }
  return 0;
}

// Action_SurfaceTension::FinishBlock()
/** Close one nblock window: S(q) = (Σ |h_q|²) / n_surfaces in the block,
  * then the same CWT γ / κ fits as Print(). Incomplete final window is
  * left in the running spectra but not written to blockout.
  */
int Action_SurfaceTension::FinishBlock() {
  if (block_surface_count_ < 1 || block_frame_count_ < 1) return 0;
  std::vector<double> spec(block_power_.size());
  for (size_t i = 0; i < spec.size(); i++)
    spec[i] = block_power_[i] / (double)block_surface_count_;
  std::vector<ST_Shell> shells;
  ST_ShellAverage(q_grid_, spec, spec, spec, shells);
  double gamma, plateau;
  double gamma_k, kappa_kT;
  int err_g = ST_CalcGamma(shells, temp_, Lt1_ref_ * Lt2_ref_, qmin_, qmax_, gamma, plateau);
  int err_k = ST_CalcKappa(shells, temp_, Lt1_ref_ * Lt2_ref_, qmin_, qmax_, gamma_k, kappa_kT);
  (void)gamma_k;
  (void)err_k;
  if (err_g) {
    mprintf("Warning: Block %i: fewer than two q shells in the fit range; skipping block gamma.\n",
            n_blocks_ + 1);
  } else {
    double wmean = block_w_sum_ / (double)block_frame_count_;
    double wtop  = block_wtop_sum_ / (double)block_frame_count_;
    double wbot  = block_wbot_sum_ / (double)block_frame_count_;
    if (block_gamma_ != 0) {
      block_gamma_->Add(n_blocks_, &gamma);
      if (block_kappa_ != 0)
        block_kappa_->Add(n_blocks_, &kappa_kT);
      block_wmean_->Add(n_blocks_, &wmean);
      block_wtop_->Add(n_blocks_, &wtop);
      block_wbot_->Add(n_blocks_, &wbot);
    }
    n_blocks_++;
  }
  std::fill(block_power_.begin(), block_power_.end(), 0.0);
  block_surface_count_ = 0;
  block_frame_count_ = 0;
  block_w_sum_ = block_wtop_sum_ = block_wbot_sum_ = 0.0;
  return 0;
}

// Action_SurfaceTension::DoAction()
/** Read the unit-cell lengths (or lx/ly/lz overrides) and ProcessFrame.
  * Skip (Action::OK) if the box is missing or a crossing failed; ERR if
  * the lateral grid changed (not NVT).
  */
Action::RetType Action_SurfaceTension::DoAction(int, ActionFrame& frm)
{
  if (!frm.Frm().BoxCrd().HasBox()) {
    SkipWarn("Frame has no box; skipping.");
    n_skipped_++;
    return Action::OK;
  }
  Vec3 lengths = frm.Frm().BoxCrd().Lengths();
  // Optional lx/ly/lz override the trajectory box (NVT slabs with noisy box records).
  double Lx = (lx_user_ > 0.0) ? lx_user_ : lengths[0];
  double Ly = (ly_user_ > 0.0) ? ly_user_ : lengths[1];
  double Lz = (lz_user_ > 0.0) ? lz_user_ : lengths[2];
  if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0) {
    SkipWarn("Invalid box lengths; skipping frame.");
    n_skipped_++;
    return Action::OK;
  }
  int err = ProcessFrame(frm.Frm(), Lx, Ly, Lz);
  if (err == 1) {
    n_skipped_++;
    return Action::OK;
  }
  if (err == 2) return Action::ERR;
  return Action::OK;
}

#ifdef MPI
// Action_SurfaceTension::SyncAction()
/** Reduce per-rank |h_q|² so Print() sees the global S(q) = ⟨|h_q|²⟩.
  * Counts (frames, surfaces, skipped) AllReduce SUM. Grid size, L₁, L₂,
  * q_min, and q_fundamental AllReduce MAX so empty ranks (zeros) do not
  * clobber. One packed ReduceMaster SUM of combined / upper / lower
  * power; |q| is a separate MAX. Roughness and block series stay on
  * DataSet::Sync (concat by rank). nblock is per-rank. Print() is master
  * only.
  */
int Action_SurfaceTension::SyncAction() {
  if (trajComm_.Size() < 2) return 0;

  int counts[4] = { n_frames_, n_surfaces_, n_skipped_, block_frame_count_ };
  trajComm_.AllReduce(counts, 4, MPI_INT, MPI_SUM);
  n_frames_          = counts[0];
  n_surfaces_        = counts[1];
  n_skipped_         = counts[2];
  block_frame_count_ = counts[3];

  int n2 = (int)total_power_.size();
  int imax[3] = { n2, nx_, ny_ };
  trajComm_.AllReduce(imax, 3, MPI_INT, MPI_MAX);
  int n2max = imax[0];
  nx_ = imax[1];
  ny_ = imax[2];
  double box[4] = { Lt1_ref_, Lt2_ref_, qmin_, q_fundamental_ };
  trajComm_.AllReduce(box, 4, MPI_DOUBLE, MPI_MAX);
  Lt1_ref_ = box[0];
  Lt2_ref_ = box[1];
  qmin_ = box[2];
  q_fundamental_ = box[3];

  if (n2max < 1) return 0;
  if (n2 != 0 && n2 != n2max) {
    rprintf("Error: surftension Fourier grid size %i differs from other ranks (%i).\n",
            n2, n2max);
    return 1;
  }
  if (n2 == 0) {
    total_power_.assign((size_t)n2max, 0.0);
    top_power_.assign((size_t)n2max, 0.0);
    bottom_power_.assign((size_t)n2max, 0.0);
    q_grid_.assign((size_t)n2max, 0.0);
  }

  // One ReduceMaster for all three spectra (latency-bound, not bandwidth).
  int npack = n2max * 3;
  std::vector<double> send((size_t)npack);
  std::copy(total_power_.begin(),  total_power_.end(),  send.begin());
  std::copy(top_power_.begin(),    top_power_.end(),    send.begin() + n2max);
  std::copy(bottom_power_.begin(), bottom_power_.end(), send.begin() + 2 * n2max);
  if (trajComm_.Master()) {
    std::vector<double> recv((size_t)npack);
    trajComm_.ReduceMaster(&recv[0], &send[0], npack, MPI_DOUBLE, MPI_SUM);
    total_power_.assign(recv.begin(), recv.begin() + n2max);
    top_power_.assign(recv.begin() + n2max, recv.begin() + 2 * n2max);
    bottom_power_.assign(recv.begin() + 2 * n2max, recv.end());
  } else {
    trajComm_.ReduceMaster(0, &send[0], npack, MPI_DOUBLE, MPI_SUM);
  }

  if (q_grid_.size() != (size_t)n2max)
    q_grid_.assign((size_t)n2max, 0.0);
  if (trajComm_.Master()) {
    std::vector<double> qmaxv((size_t)n2max);
    trajComm_.ReduceMaster(&qmaxv[0], &q_grid_[0], n2max, MPI_DOUBLE, MPI_MAX);
    q_grid_.swap(qmaxv);
    grid_ready_ = true;
  } else {
    trajComm_.ReduceMaster(0, &q_grid_[0], n2max, MPI_DOUBLE, MPI_MAX);
  }
  return 0;
}
#endif

// Action_SurfaceTension::Print()
/** Form S(q) = (Σ |h_q|²) / n_surfaces, shell-average, then the CWT fits
  *   γ = k_B T / (A ⟨q² S⟩) ,   1/(q² S) = (A / k_B T) (γ + κ q²)
  * Fill the q-meshes, print the numeric summary, and write summaryout.
  * Per-block γ / κ are printed here, not during the frame loop.
  */
void Action_SurfaceTension::Print()
{
  const char* t1 = 0;
  const char* t2 = 0;
  ST_PlaneAxes((int)normal_, t1, t2);
  mprintf("    SURFTENSION:\n");
  mprintf("\tSlab normal: %s (interface plane %s-%s)\n",
          ST_AxisName((int)normal_), t1, t2);
  if (iface_ == WILLARD)
    mprintf("\tInterface: Willard-Chandler density isosurface.\n");
  else
    mprintf("\tInterface: ITIM min/max.\n");
  mprintf("\tSurfaces: %i", nsurf_);
  if (nsurf_ == 1)
    mprintf(" (%s)", (side_ == SIDE_UPPER) ? "upper" : "lower");
  mprintf("\n");
  if (has_mask2_)
    mprintf("\tMask2 (lower): '%s'\n", Mask2_.MaskString());
  if (nblock_ > 0 && block_frame_count_ > 0) {
    mprintf("\tNOTE: Incomplete nblock window (%i frames) excluded from blockout.\n",
            block_frame_count_);
  }
  if (n_frames_ < 1) {
    mprinterr("Error: surftension: No frames were analyzed");
    if (n_skipped_ > 0)
      mprinterr(" (%i skipped)", n_skipped_);
    mprinterr(".\n");
    if (summaryFile_ != 0) {
      summaryFile_->Printf("# surftension summary\n");
      ST_SummaryI(summaryFile_, "frames", 0);
      ST_SummaryI(summaryFile_, "skipped", n_skipped_);
    }
    return;
  }

  std::vector<double> spec(total_power_.size());
  std::vector<double> spec_top(top_power_.size());
  std::vector<double> spec_bot(bottom_power_.size());
  for (size_t i = 0; i < spec.size(); i++) {
    spec[i] = total_power_[i] / (double)n_surfaces_;
    spec_top[i] = top_power_[i] / (double)n_frames_;
    spec_bot[i] = bottom_power_[i] / (double)n_frames_;
  }
  std::vector<ST_Shell> shells;
  ST_ShellAverage(q_grid_, spec, spec_top, spec_bot, shells);

  // Fit γ on the combined, upper-only, and lower-only shells.
  double area = Lt1_ref_ * Lt2_ref_;
  double gamma_full, plateau;
  double gamma_top, plateau_top;
  double gamma_bot, plateau_bot;
  std::vector<ST_Shell> top_only = shells;
  std::vector<ST_Shell> bot_only = shells;
  for (size_t i = 0; i < shells.size(); i++) {
    top_only[i].S = shells[i].Stop;
    bot_only[i].S = shells[i].Sbot;
  }
  int err_full = ST_CalcGamma(shells, temp_, area, qmin_, qmax_, gamma_full, plateau);
  int err_top = 1, err_bot = 1;
  if (do_upper_)
    err_top = ST_CalcGamma(top_only, temp_, area, qmin_, qmax_, gamma_top, plateau_top);
  if (do_lower_)
    err_bot = ST_CalcGamma(bot_only, temp_, area, qmin_, qmax_, gamma_bot, plateau_bot);
  (void)plateau_top;
  (void)plateau_bot;
  double gamma_h, kappa_h, gamma_h_top, kappa_h_top, gamma_h_bot, kappa_h_bot;
  int err_kh     = ST_CalcKappa(shells, temp_, area, qmin_, qmax_, gamma_h, kappa_h);
  int err_kh_top = 1, err_kh_bot = 1;
  if (do_upper_)
    err_kh_top = ST_CalcKappa(top_only, temp_, area, qmin_, qmax_, gamma_h_top, kappa_h_top);
  if (do_lower_)
    err_kh_bot = ST_CalcKappa(bot_only, temp_, area, qmin_, qmax_, gamma_h_bot, kappa_h_bot);
  double gamma_for_kappaq = err_kh ? gamma_full : gamma_h;

  double slope     = ST_LogSlope(shells, qmin_, qmax_, &ST_Shell::S);
  double slope_top = ST_LogSlope(shells, qmin_, qmax_, &ST_Shell::Stop);
  double slope_bot = ST_LogSlope(shells, qmin_, qmax_, &ST_Shell::Sbot);

  for (size_t i = 0; i < shells.size(); i++) {
    double q = shells[i].q;
    S_->AddXY(q, shells[i].S);
    S_top_->AddXY(q, shells[i].Stop);
    S_bot_->AddXY(q, shells[i].Sbot);
    q2S_->AddXY(q, q * q * shells[i].S);
    q2S_top_->AddXY(q, q * q * shells[i].Stop);
    q2S_bot_->AddXY(q, q * q * shells[i].Sbot);
    gammaq_->AddXY(q, ST_ShellGamma(q, shells[i].S, temp_, area));
    gammaq_top_->AddXY(q, ST_ShellGamma(q, shells[i].Stop, temp_, area));
    gammaq_bot_->AddXY(q, ST_ShellGamma(q, shells[i].Sbot, temp_, area));
    kappaq_->AddXY(q, ST_ShellKappa(q, shells[i].S, temp_, area, gamma_for_kappaq));
    kappaq_top_->AddXY(q, ST_ShellKappa(q, shells[i].Stop, temp_, area,
                                       err_kh_top ? gamma_top : gamma_h_top));
    kappaq_bot_->AddXY(q, ST_ShellKappa(q, shells[i].Sbot, temp_, area,
                                       err_kh_bot ? gamma_bot : gamma_h_bot));
  }

  double mean_w = 0.0, mean_wt = 0.0, mean_wb = 0.0;
  if (wmean_->Size() > 0) {
    int nu = 0, nl = 0;
    for (size_t i = 0; i < wmean_->Size(); i++) {
      mean_w += ((DataSet_1D*)wmean_)->Dval(i);
      if (do_upper_) {
        mean_wt += ((DataSet_1D*)wtop_)->Dval(i);
        nu++;
      }
      if (do_lower_) {
        mean_wb += ((DataSet_1D*)wbot_)->Dval(i);
        nl++;
      }
    }
    mean_w /= (double)wmean_->Size();
    if (nu > 0) mean_wt /= (double)nu;
    if (nl > 0) mean_wb /= (double)nl;
  }

  mprintf("\tT = %g K\n", temp_);
  mprintf("\tFrames analyzed = %i", n_frames_);
  if (n_skipped_ > 0)
    mprintf(" (%i skipped)", n_skipped_);
  mprintf("\n");
  mprintf("\tL(%s) = %g Ang, L(%s) = %g Ang (from unit cell unless lx/ly/lz set)\n",
          t1, Lt1_ref_, t2, Lt2_ref_);
  mprintf("\tArea = %g Ang^2\n", area);
  if (q_fundamental_ > 0.0)
    mprintf("\tFundamental |q| = %g Ang^-1\n", q_fundamental_);
  mprintf("\tFit q range = %g to %g Ang^-1\n", qmin_, qmax_);
  if (ST_Finite(slope))
    mprintf("\tLow-q log-log slope = %g (ideal capillary-wave slope is about -2)\n", slope);
  if (err_full)
    mprinterr("Error: Combined-surface gamma fit failed.\n");
  else
    mprintf("\tgamma, combined surfaces = %g mN/m (q^2 S plateau = %g)\n", gamma_full, plateau);
  if (!err_top)
    mprintf("\tgamma, upper surface = %g mN/m\n", gamma_top);
  if (!err_bot)
    mprintf("\tgamma, lower surface = %g mN/m\n", gamma_bot);
  if (err_kh)
    mprintf("\tNOTE: Helfrich kappa fit needs at least three q shells in the fit range.\n");
  else {
    mprintf("\tgamma, Helfrich intercept = %g mN/m\n", gamma_h);
    mprintf("\tkappa, combined surfaces = %g kT (%g J)\n",
            kappa_h, kappa_h * ST_KB * temp_);
  }
  if (!err_kh_top)
    mprintf("\tkappa, upper surface = %g kT\n", kappa_h_top);
  if (!err_kh_bot)
    mprintf("\tkappa, lower surface = %g kT\n", kappa_h_bot);
  if (do_upper_ && do_lower_ && ST_Finite(slope_top) && ST_Finite(slope_bot))
    mprintf("\tLow-q slope upper/lower = %g / %g\n", slope_top, slope_bot);
  if (do_upper_ && do_lower_)
    mprintf("\tMean roughness = %g Ang (upper %g, lower %g)\n", mean_w, mean_wt, mean_wb);
  else if (do_upper_)
    mprintf("\tMean roughness (upper) = %g Ang\n", mean_wt);
  else
    mprintf("\tMean roughness (lower) = %g Ang\n", mean_wb);

  if (block_gamma_ != 0 && block_gamma_->Size() > 0) {
    mprintf("\tCompleted blocks = %zu (nblock = %i frames)\n",
            block_gamma_->Size(), nblock_);
    for (size_t i = 0; i < block_gamma_->Size(); i++) {
      double g = ((DataSet_1D*)block_gamma_)->Dval(i);
      double w = (block_wmean_ != 0) ? ((DataSet_1D*)block_wmean_)->Dval(i) : ST_NaN();
      double k = ST_NaN();
      if (block_kappa_ != 0 && i < block_kappa_->Size())
        k = ((DataSet_1D*)block_kappa_)->Dval(i);
      if (ST_Finite(k))
        mprintf("\tBlock %zu: gamma = %g mN/m, kappa = %g kT, roughness = %g Ang\n",
                i + 1, g, k, w);
      else
        mprintf("\tBlock %zu: gamma = %g mN/m, roughness = %g Ang\n",
                i + 1, g, w);
    }
  }
  double block_gmean = ST_NaN(), block_gsd = ST_NaN(), block_gsem = ST_NaN();
  double block_kmean = ST_NaN(), block_ksd = ST_NaN(), block_ksem = ST_NaN();
  if (block_gamma_ != 0 && block_gamma_->Size() > 1) {
    double gmean = 0.0, g2 = 0.0;
    for (size_t i = 0; i < block_gamma_->Size(); i++) {
      double g = ((DataSet_1D*)block_gamma_)->Dval(i);
      gmean += g;
      g2 += g * g;
    }
    gmean /= (double)block_gamma_->Size();
    double var = (g2 - (double)block_gamma_->Size() * gmean * gmean) /
                 (double)(block_gamma_->Size() - 1);
    double sd = (var > 0.0) ? sqrt(var) : 0.0;
    double sem = sd / sqrt((double)block_gamma_->Size());
    block_gmean = gmean;
    block_gsd = sd;
    block_gsem = sem;
    mprintf("\tBlock mean gamma = %g mN/m\n", gmean);
    mprintf("\tBlock SD gamma = %g mN/m\n", sd);
    mprintf("\tBlock SEM gamma = %g mN/m\n", sem);
    mprintf("\tgamma +/- 2 SEM = %g mN/m\n", 2.0 * sem);
  }
  if (block_kappa_ != 0 && block_kappa_->Size() > 1) {
    double kmean = 0.0, k2 = 0.0;
    int nk = 0;
    for (size_t i = 0; i < block_kappa_->Size(); i++) {
      double k = ((DataSet_1D*)block_kappa_)->Dval(i);
      if (!ST_Finite(k)) continue;
      kmean += k;
      k2 += k * k;
      nk++;
    }
    if (nk > 1) {
      kmean /= (double)nk;
      double var = (k2 - (double)nk * kmean * kmean) / (double)(nk - 1);
      double sd = (var > 0.0) ? sqrt(var) : 0.0;
      double sem = sd / sqrt((double)nk);
      block_kmean = kmean;
      block_ksd = sd;
      block_ksem = sem;
      mprintf("\tBlock mean kappa = %g kT\n", kmean);
      mprintf("\tBlock SD kappa = %g kT\n", sd);
      mprintf("\tBlock SEM kappa = %g kT\n", sem);
      mprintf("\tkappa +/- 2 SEM = %g kT\n", 2.0 * sem);
    }
  }

  if (summaryFile_ != 0) {
    summaryFile_->Printf("# surftension summary\n");
    ST_SummaryD(summaryFile_, "temperature", temp_);
    ST_SummaryI(summaryFile_, "nsurf", nsurf_);
    ST_SummaryS(summaryFile_, "normal", ST_AxisName((int)normal_));
    ST_SummaryS(summaryFile_, "interface", (iface_ == WILLARD) ? "willard" : "itim");
    ST_SummaryS(summaryFile_, "mask", Mask_.MaskString());
    if (has_mask2_)
      ST_SummaryS(summaryFile_, "mask2", Mask2_.MaskString());
    ST_SummaryD(summaryFile_, "L_t1", Lt1_ref_);
    ST_SummaryD(summaryFile_, "L_t2", Lt2_ref_);
    ST_SummaryD(summaryFile_, "area", area);
    ST_SummaryD(summaryFile_, "q_fundamental", q_fundamental_);
    ST_SummaryD(summaryFile_, "qmin", qmin_);
    ST_SummaryD(summaryFile_, "qmax", qmax_);
    ST_SummaryI(summaryFile_, "frames", n_frames_);
    ST_SummaryI(summaryFile_, "skipped", n_skipped_);
    if (!err_full)
      ST_SummaryD(summaryFile_, "gamma", gamma_full);
    if (do_upper_ && !err_top)
      ST_SummaryD(summaryFile_, "gamma_upper", gamma_top);
    if (do_lower_ && !err_bot)
      ST_SummaryD(summaryFile_, "gamma_lower", gamma_bot);
    if (!err_kh)
      ST_SummaryD(summaryFile_, "kappa", kappa_h);
    ST_SummaryD(summaryFile_, "roughness", mean_w);
    if (do_upper_)
      ST_SummaryD(summaryFile_, "roughness_upper", mean_wt);
    if (do_lower_)
      ST_SummaryD(summaryFile_, "roughness_lower", mean_wb);
    ST_SummaryD(summaryFile_, "block_mean_gamma", block_gmean);
    ST_SummaryD(summaryFile_, "block_sd_gamma", block_gsd);
    ST_SummaryD(summaryFile_, "block_sem_gamma", block_gsem);
    ST_SummaryD(summaryFile_, "block_mean_kappa", block_kmean);
    ST_SummaryD(summaryFile_, "block_sd_kappa", block_ksd);
    ST_SummaryD(summaryFile_, "block_sem_kappa", block_ksem);
  }
}
