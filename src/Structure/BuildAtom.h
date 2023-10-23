#ifndef INC_STRUCTURE_BUILDATOM_H
#define INC_STRUCTURE_BUILDATOM_H
#include "StructureEnum.h"
namespace Cpptraj {
namespace Structure {
/// Hold information for an atom used when building/modelling new coordinates.
class BuildAtom {
  public:
    BuildAtom() : positionKnown_(false), ctype_(IS_UNKNOWN_CHIRALITY) {}
    /// Set position status
    void SetPositionKnown(bool b) { positionKnown_ = b; }
    /// Set chirality
    void SetChirality(ChiralType c) { ctype_ = c; }
    /// Set size of priority array based on number of bonds
    void SetNbonds(int n) { priority_.assign(n, -1); }
    /// \return Pointer to priority array
    int* PriorityPtr() { return &(priority_[0]); }
  private:
    bool positionKnown_;        ///< True if position is "known", i.e. can be used to build.
    ChiralType ctype_;          ///< Chirality around atom
    std::vector<int> priority_; ///< Indices of bonded atoms, sorted by priority
};
}
}
#endif
