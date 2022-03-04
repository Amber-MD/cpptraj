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
    SetFromFrame(Frame const&, Topology const&);
  private:
    typedef std::vector<InternalCoords> ICarray;

    ICarray IC_; ///< Hold internal coordinates for all atoms
};
}
}
#endif
