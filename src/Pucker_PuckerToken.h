#ifndef INC_PUCKER_PUCKERTOKEN_H
#define INC_PUCKER_PUCKERTOKEN_H
#include <vector>
#include "NameType.h"
namespace Cpptraj {
namespace Pucker {
/// Used to search for atoms involved in puckering in a Topology
class PuckerToken {
  public:
    PuckerToken();

    typedef std::vector<NameType> NameArray;
  private:
    NameArray atomNames_; ///< Atoms that define the pucker.
};

}
}
#endif
