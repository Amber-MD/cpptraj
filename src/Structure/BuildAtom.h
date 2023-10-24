#ifndef INC_STRUCTURE_BUILDATOM_H
#define INC_STRUCTURE_BUILDATOM_H
#include "StructureEnum.h"
namespace Cpptraj {
namespace Structure {
/// Hold information for an atom used when building/modelling new coordinates.
class BuildAtom {
  public:
    BuildAtom() : ctype_(IS_UNKNOWN_CHIRALITY) {}
    /// Set chirality
    void SetChirality(ChiralType c) { ctype_ = c; }
    /// Set size of priority array based on number of bonds
    //void SetNbonds(int n) { priority_.assign(n, -1); }
    /// \return Pointer to priority array
    //int* PriorityPtr() { return &(priority_[0]); }
    /// Used to modify the priority array
    std::vector<int>& ModifyPriority() { return priority_; }

    /// \return Atom chirality
    ChiralType Chirality() const { return ctype_; }
    /// \return Priority array
    std::vector<int> const& Priority() const { return priority_; }
  private:
    ChiralType ctype_;          ///< Chirality around atom
    std::vector<int> priority_; ///< Indices of bonded atoms, sorted by priority
};
}
}
#endif
