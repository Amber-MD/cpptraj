#ifndef INC_ENERGY_CMAP_H
#define INC_ENERGY_CMAP_H
#include <vector>
#include "../Matrix.h"
class CharMask;
class CmapArray;
class CmapGridArray;
class Frame;
class Topology;
class Vec3;
namespace Cpptraj {
namespace Energy {
/// Implement CMAP energy term
class CMAP {
  public:
    CMAP();
    ~CMAP();

    /// Setup CMAP splines etc
    int Setup_CMAP_Ene(Topology const&, CharMask const&);
    /// Calculate CMAP energy
    double Ene_CMAP(Frame const&) const;
    /// Calculate CMAP energy and force
    double Ene_Frc_CMAP(Frame const&) const;
  private:
    //double get_cmap_energy(CmapArray const&, Frame const&, double&, double&,
    //                         Vec3(&dPhi_dijkl)[4], Vec3(&dPsi_djklm)[4]) const;
    /// Calculate CMAP energy
    double Ene_CMAP(CmapArray const&, Frame const&) const;
    /// Calculate CMAP energy and force
    double Ene_Frc_CMAP(CmapArray const&, Frame const&) const;

    double charmm_calc_cmap_from_phi_psi(double, double, int, double&, double&) const;
    static double evaluate_cubic_spline(int, std::vector<double> const&, std::vector<double> const&, int);
    int generate_cmap_derivatives(Topology const&);
    static void generate_cubic_spline(int, int, std::vector<double> const&, std::vector<double>&);

    typedef Matrix<double> Dmatrix;
    typedef std::vector<Dmatrix> DMarray;

    DMarray cmap_dPhi_; ///< Hold dE/dPhi for each CMAP grid
    DMarray cmap_dPsi_; ///< Hold dE/dPsi for each CMAP grid
    DMarray cmap_dPhi_dPsi_; ///< Hold d^2E/dPhidPsi for each CMAP grid
    CmapGridArray const* cmapGridPtr_; ///< Pointer to CMAP grids this was set up with.
    CmapArray* selected_cmaps_; ///< Hold cmaps to calculate energy for if not all selected.
    CmapArray const* all_cmaps_; ///< Pointer to CMAP array in Topology when all selected.

    static const int wt_[16][16]; ///< Weight matrix
};
}
}
#endif
