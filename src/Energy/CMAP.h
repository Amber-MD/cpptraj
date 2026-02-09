#ifndef INC_ENERGY_CMAP_H
#define INC_ENERGY_CMAP_H
#include <vector>
class CmapArray;
class Frame;
class Topology;
namespace Cpptraj {
namespace Energy {
/// Implement CMAP energy term
class CMAP {
  public:
    CMAP();

    /// Calculate CMAP energy
    double Ene_CMAP(CmapArray const&, Frame const&) const;
  private:
    static double evaluate_cubic_spline(int, std::vector<double> const&, std::vector<double> const&, int);
    int generate_cmap_derivatives(Topology const&);
};
}
}
#endif
