#ifndef INC_ENERGY_CMAP_H
#define INC_ENERGY_CMAP_H
class CmapArray;
class Frame;
namespace Cpptraj {
namespace Energy {
/// Implement CMAP energy term
class CMAP {
  public:
    CMAP();

    /// Calculate CMAP energy
    double Ene_CMAP(CmapArray const&, Frame const&) const;
  private:
};
}
}
#endif
