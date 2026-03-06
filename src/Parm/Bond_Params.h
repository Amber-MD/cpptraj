#ifndef INC_PARM_BOND_PARAMS_H
#define INC_PARM_BOND_PARAMS_H
#include "../Atom.h"
#include <vector>
#include <set>
class Topology;
class BondType;
class BondParmArray;
namespace Cpptraj {
namespace Parm {
typedef std::vector< std::set<Atom::AtomicElementType> > BP_mapType;
/// Generate bond parameters
void GenerateBondParam(BondParmArray&, BondType&, BP_mapType&, std::vector<Atom> const&);
/// Used to generate bond parameters (Req only) based on atomic element types.
void GenerateBondParameters(Topology&);
}
}
#endif
