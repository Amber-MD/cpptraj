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
  private:
    typedef std::vector<Cpptraj::Chirality::ChiralType> Carray;
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> Parray;

    Carray atomChirality_; ///< Hold chirality for each atom in the combined system.
    Parray atomPriority_;  ///< Hold bonded atom index priority array for each atom.
};
}
}
#endif
