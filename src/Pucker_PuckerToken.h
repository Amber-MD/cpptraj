#ifndef INC_PUCKER_PUCKERTOKEN_H
#define INC_PUCKER_PUCKERTOKEN_H
#include <vector>
#include "NameType.h"
// Forward declares
class Topology;
namespace Cpptraj {
namespace Pucker {
class PuckerMask;
/// Used to define atoms involved in puckering
class PuckerToken {
  public:
    PuckerToken();
    typedef std::vector<NameType> NameArray;
    /// CONSTRUCTOR - take name and array of atom names
    PuckerToken(std::string const&, NameArray const&);
    /// \return pucker token name
    std::string const& Name() const { return name_; }

    PuckerMask FindPuckerAtoms(Topology const&, int) const;
  private:
    void FindAtoms(Topology const&, int, unsigned int, unsigned int, std::vector<int>&) const;

    std::string name_;    ///< Pucker name.
    NameArray atomNames_; ///< Atoms that define the pucker.
};

}
}
#endif
