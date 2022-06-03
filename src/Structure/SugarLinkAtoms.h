#ifndef INC_STRUCTURE_SUGARLINKATOMS_H
#define INC_STRUCTURE_SUGARLINKATOMS_H
#include "LinkAtom.h"
#include <set>
#include <string>
// Forward declares
class Topology;
namespace Cpptraj {
namespace Structure {
/// Hold sugar link atoms for determining overall linkage.
class SugarLinkAtoms {
  public:
    /// CONSTRUCTOR - takes debug level
    SugarLinkAtoms(int);
    /// Add link atom with given topology index and chain position
    void AddLinkAtom(int at, int pos) { linkages_.insert(LinkAtom(at, pos)); }
    /// \return True if no link atoms present
    bool NoLinkAtoms() const { return linkages_.empty(); }
    /// \return Glycam linkage code for given linked atoms
    std::string GlycamLinkageCode(Topology const&) const;
  private:
    std::set<LinkAtom> linkages_;
    int debug_;
};

}
}
#endif
