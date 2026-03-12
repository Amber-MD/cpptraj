#ifndef INC_STRUCTURE_OLDBUILDER_H
#define INC_STRUCTURE_OLDBUILDER_H
#include <vector>
#include "StructureEnum.h"
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
/// This is the old (<v7) builder used to attach different topology/frame combos using internal coordinates.
/** Retaining this class for backwards compatibility with modXNA. */
class OldBuilder {
  public:
    /// CONSTRUCTOR
    OldBuilder();
    /// Combine second fragment into first fragment and bond
    int Combine(Topology&, Frame&, Topology const&, Frame const&, int, int);
    /// Set debug
    void SetDebug(int d) { debug_ = d; }
  private:
    typedef std::vector<bool> Barray;

    int debug_;
};

/// Old (<v7) class for holding information for an atom used when building/modelling new coordinates.
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
