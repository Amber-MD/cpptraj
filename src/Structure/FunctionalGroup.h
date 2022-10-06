#ifndef INC_STRUCTURE_FUNCTIONALGROUP_H
#define INC_STRUCTURE_FUNCTIONALGROUP_H
#include "../Atom.h"
#include "../NameType.h"
#include <vector>
class Topology;
namespace Cpptraj {
namespace Structure {
/// Hold information for a functional group TODO use this to store functional group data
class FunctionalGroup {
  public:
    /// Different functional group types
    enum Type {
      G_SO3 = 0, ///< Sulfate
      G_CH3,     ///< Methyl
      G_ACX,     ///< Acetyl
      G_OH,      ///< Hydroxyl
      G_OME,     ///< O-methyl
      UNRECOGNIZED_GROUP };
    /// CONSTRUCTOR TODO remove blank constructor?
    FunctionalGroup();
    /// CONSTRUCTOR - take name
    //FunctionalGroup(NameType const&);

    /// Set up group from given topology.
    int SetupFromTop(Topology const&, Atom::AtomicElementType);

    /// Print info to stdout
    void PrintInfo() const;

    /// \return true if given FunctionalGroup atom IDs match this ones
    bool Match(FunctionalGroup const&) const;

    static const char* typeString(Type t) { return FunctionalGroupStr_[t]; }
  private:
    /// Keep synced with FunctionalGroupType
    static const char* FunctionalGroupStr_[];

    /// Clear all information
    void Clear();

    /// Add atom and ID to group
    void AddAtom(NameType const&, std::string const&);

    NameType resname_;                    ///< Functional group residue name.
    std::vector<NameType> anames_;        ///< Functional group atom names. Heavy atoms first.
    std::vector<std::string> atomIDs_;    ///< Functional group atom IDs.
    Atom::AtomicElementType linkAtomElt_; ///< Element of atom group is linked to.
    Atom::AtomicElementType chargeAtom_;  ///< Element of atom which needs charge adjusted.
    double chargeOffset_;                 ///< Charge offset for adjusting charge.
};

}
}
#endif
