#ifndef INC_STRUCTURE_CHIRALITY_H
#define INC_STRUCTURE_CHIRALITY_H
#include <vector>
#include "StructureEnum.h"
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {

/// \return Chirality at specified atom, set torsion value
ChiralType DetermineChirality(double&, int*, int, Topology const&, Frame const&, int);
/// \return Chirality at specified atom
ChiralType DetermineChirality(int, Topology const&, Frame const&, int);
/// \return Chirality at specified atom, set priority
ChiralType SetPriority(std::vector<int>&, int, Topology const&, Frame const&, int);
}
}
#endif
