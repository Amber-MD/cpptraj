#ifndef INC_STRUCTURE_ZMATRIX_H
#define INC_STRUCTURE_ZMATRIX_H
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
class InternalCoords;
/// Hold internal coordinates for a system.
class Zmatrix {
  public:
    /// CONSTRUCTOR
    Zmatrix();
    /// Convert Frame/Topology to internal coordinates array
    int SetFromFrame(Frame const&, Topology const&);
  private:
    typedef std::vector<InternalCoords> ICarray;

    ICarray IC_; ///< Hold internal coordinates for all atoms
    int seed0_;  ///< Index of first seed atom
    int seed1_;  ///< Index of second seed atom
    int seed2_;  ///< Index of third seed atom
};
}
}
#endif
