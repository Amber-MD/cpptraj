#ifndef INC_ENERGY_CMAP_H
#define INC_ENERGY_CMAP_H
#include <vector>
#include "../Matrix.h"
class CmapArray;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
/// Implement CMAP energy term
class CMAP {
  public:
    CMAP();

    /// Setup CMAP splines etc
    int Setup_CMAP_Ene(Topology const&);
    /// Calculate CMAP energy
    double Ene_CMAP(CmapArray const&, Frame const&) const;
  private:
    static double evaluate_cubic_spline(int, std::vector<double> const&, std::vector<double> const&, int);
    int generate_cmap_derivatives(Topology const&);
    static void generate_cubic_spline(int, int, std::vector<double> const&, std::vector<double>&);

    typedef Matrix<double> Dmatrix;
    typedef std::vector<Dmatrix> DMarray;

    DMarray cmap_dPhi_; ///< Hold dE/dPhi for each CMAP grid
    DMarray cmap_dPsi_; ///< Hold dE/dPsi for each CMAP grid
    DMarray cmap_dPhi_dPsi_; ///< Hold d^2E/dPhidPsi for each CMAP grid
};
}
}
#endif
