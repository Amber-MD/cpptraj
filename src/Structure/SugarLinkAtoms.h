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

    /// Iterator over link atoms
    typedef std::set<LinkAtom>::const_iterator const_iterator;
    /// \return iterator to beginning of link atoms
    const_iterator begin() const { return linkages_.begin(); }
    /// \return iterator to end of link atoms
    const_iterator end()   const { return linkages_.end(); }
  private:
    std::set<LinkAtom> linkages_;
    int debug_;
};

}
}
#endif
