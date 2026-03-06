#ifndef INC_STRUCTURE_CHIRALITY_H
#define INC_STRUCTURE_CHIRALITY_H
#include <vector>
#include "StructureEnum.h"
class Atom;
class Frame;
class Topology;
class Vec3;
namespace Cpptraj {
namespace Structure {

/// \return Chirality at specified atom, set torsion value
ChiralType DetermineChirality(double&, int*, int, Topology const&, Frame const&, int);
/// \return Chirality at specified atom
ChiralType DetermineChirality(int, Topology const&, Frame const&, int);
/// \return Chirality at specified atom, set priority
ChiralType SetPriority(std::vector<int>&, int, Topology const&, Frame const&, int);
/// \return Chirality at specified atom, set priority, no warnings
ChiralType SetPriority_silent(std::vector<int>&, int, Topology const&, Frame const&);

/// This namespace will hold the LEaP chirality routines
namespace Chirality {

/// Calculate chirality in same manner as LEaP. All atom positions should be known.
double VectorAtomChirality(Vec3 const&, Vec3 const&, Vec3 const&, Vec3 const&);
/// Calculate chirality in same manner as LEaP. Some atom positions may not be known.
double VectorAtomNormalizedChirality(Vec3 const&,
                                     Vec3 const&, bool, Vec3 const&, bool,
                                     Vec3 const&, bool, Vec3 const&, bool);
/// Order atoms for chirality calculation like LEaP
void chiralityOrderNeighbors(Atom const&, int&, int&, int&, int&);
/// Transform given chirality to an orientation
double chiralityToOrientation(double, Atom const&, int, int, int, int);

}
}
}
#endif
