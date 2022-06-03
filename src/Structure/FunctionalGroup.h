#ifndef INC_STRUCTURE_FUNCTIONALGROUP_H
#define INC_STRUCTURE_FUNCTIONALGROUP_H
#include "../Atom.h"
#include "../NameType.h"
#include <vector>
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

    FunctionalGroup();

    static const char* typeString(Type t) { return FunctionalGroupStr_[t]; }
  private:
    /// Keep synced with FunctionalGroupType
    static const char* FunctionalGroupStr_[];

    NameType resname_;                   ///< Functional group residue name.
    std::vector<NameType> anames_;       ///< Functional group atom names. Heavy atoms first.
    Atom::AtomicElementType chargeAtom_; ///< Element of atom which needs charge adjusted.
    double chargeOffset_;                ///< Charge offset for adjusting charge.
};

}
}
#endif
