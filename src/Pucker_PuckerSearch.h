#ifndef INC_PUCKER_PUCKERSEARCH_H
#define INC_PUCKER_PUCKERSEARCH_H
#include <vector>
#include "Pucker_PuckerToken.h"
#include "Pucker_PuckerMask.h"
class Range;
class Topology;
namespace Cpptraj {
namespace Pucker {
// Forward declares
class PuckerToken;
/// Used to search for puckers in a Topology
class PuckerSearch {
  public:
    PuckerSearch();
    /// Search for any defined puckers in residue range in given Topology
    int FindPuckers(Topology const&, Range const&);
  private:
    std::vector<PuckerToken> puckersToSearchFor_; ///< List of puckers to search for
    std::vector<PuckerMask> foundPuckers_;        ///< List of found puckers
};

}
}
#endif
