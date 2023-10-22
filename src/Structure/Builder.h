#ifndef INC_STRUCTURE_BUILDER_H
#define INC_STRUCTURE_BUILDER_H
#include <vector>
#include "../Chirality.h"
namespace Cpptraj {
namespace Structure {
/// Used to attach different topology/frame combos using internal coordinates
class Builder {
  public:
    /// CONSTRUCTOR
    Builder();
    /// Combine smaller fragment into larger fragment
    int Combine(Topology&, Frame&, Topology&, Frame&, int, int);
  private:
    typedef std::vector<Cpptraj::Chirality::ChiralType> Carray;
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> Parray;
    typedef std::vector<bool> Barray;

    /// \return heavy atom count
    static inline int heavy_atom_count(Topology const&);
    /// Combine fragment1 into fragment 0
    int combine01(Topology&, Frame&, Topology const&, Frame const&, int, int);

    Carray atomChirality_; ///< Hold chirality for each atom in the combined system.
    Parray atomPriority_;  ///< Hold bonded atom index priority array for each atom.
    int debug_;
};
}
}
#endif
