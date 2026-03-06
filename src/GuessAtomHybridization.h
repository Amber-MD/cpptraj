#ifndef INC_GUESSATOMHYBRIDIZATION_H
#define INC_GUESSATOMHYBRIDIZATION_H
#include <vector>
#include "AtomType.h" // for HybridizationType
class Atom;
namespace Cpptraj {
/// Guess atom hybridization type based on element/bonding
AtomType::HybridizationType GuessAtomHybridization(Atom const&, std::vector<Atom> const&);
}
#endif
