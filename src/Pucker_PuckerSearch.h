#ifndef INC_PUCKER_PUCKERSEARCH_H
#define INC_PUCKER_PUCKERSEARCH_H
#include <vector>
#include "Pucker_PuckerToken.h"
#include "Pucker_PuckerMask.h"
#include "Pucker.h"
class Range;
class Topology;
class ArgList;
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
    /// Indicate we want to search for the given pre-defined type
    int SearchFor(Type);
    /// Look at ArgList for recognized pucker keywords
    int SearchForArgs(ArgList&);
    /// Search for new types defined in ArgList
    int SearchForNewTypeArgs(ArgList&);
    /// Search for all types if no types yet defined
    int SearchForAll();
  private:
    static const char* Keywords_[];               ///< Keywords corresponding to Pucker::Type
    /// Define a new custom pucker type
    int SearchForNewType(std::string const&, PuckerToken::NameArray const&);

    std::vector<PuckerToken> puckersToSearchFor_; ///< List of puckers to search for
    std::vector<PuckerMask> foundPuckers_;        ///< List of found puckers
};

}
}
#endif
