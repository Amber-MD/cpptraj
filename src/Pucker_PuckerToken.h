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
    PuckerToken(NameArray const&);

    PuckerMask FindPuckerAtoms(Topology const&, int) const;
  private:
    void FindAtoms(Topology const&, int, unsigned int, unsigned int, std::vector<int>&) const;

    NameArray atomNames_; ///< Atoms that define the pucker.
};

}
}
#endif
