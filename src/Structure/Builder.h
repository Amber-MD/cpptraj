#ifndef INC_STRUCTURE_BUILDER_H
#define INC_STRUCTURE_BUILDER_H
#include <vector>
#include "BuildAtom.h"
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
class BuildAtom;
class Zmatrix;
/// Used to attach different topology/frame combos using internal coordinates
class Builder {
  public:
    /// CONSTRUCTOR
    Builder();
    /// Combine smaller fragment into larger fragment
    int Combine(Topology&, Frame&, Topology const&, Frame const&, Topology const&, Frame const&, int, int);
  private:
    typedef std::vector<BuildAtom> AtArray;
    typedef std::vector<bool> Barray;

    /// \return heavy atom count
    static inline int heavy_atom_count(Topology const&);
    /// Combine fragment1 into fragment 0
    int combine01(Topology&, Frame&, Topology const&, Frame const&, Topology const&, Frame const&, int, int);

    AtArray atoms_; ///< Hold chirality, bonded atom index priority, and position status for each atom
    int debug_;
};
}
}
#endif
