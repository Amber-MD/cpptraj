#ifndef INC_STRUCTURE_ZMATRIX_H
#define INC_STRUCTURE_ZMATRIX_H
#include "InternalCoords.h"
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
/// Hold internal coordinates for a system.
class Zmatrix {
    typedef std::vector<InternalCoords> ICarray;
  public:
    /// CONSTRUCTOR
    Zmatrix();
    /// Add internal coordinate
    void AddIC(InternalCoords const&);
    /// Add internal coordinate as next available seed
    int AddICseed(InternalCoords const&);
    /// Convert Frame/Topology to internal coordinates array
    int SetFromFrame(Frame const&, Topology const&);
    /// Set Frame from internal coords
    int SetToFrame(Frame&) const;
    /// Print to stdout
    void print() const;

    typedef ICarray::const_iterator const_iterator;
    const_iterator begin() const { return IC_.begin(); }
    const_iterator end()   const { return IC_.end(); }
  private:

    ICarray IC_; ///< Hold internal coordinates for all atoms
    int seed0_;  ///< Index of first seed atom
    int seed1_;  ///< Index of second seed atom
    int seed2_;  ///< Index of third seed atom
};
}
}
#endif
