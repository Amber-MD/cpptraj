#ifndef INC_STRUCTURE_MODEL_H
#define INC_STRUCTURE_MODEL_H
#include <vector>
#include "StructureEnum.h"
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
class BuildAtom;
/// Routines to generate model parameters
namespace Model {
/// Given atoms J K and L, attempt to assign a reasonable value for phi for atom I
int AssignPhi(double&, int, int, int, int, Topology const&, Frame const&, std::vector<bool> const&, BuildAtom const&, int);
/// Given atoms J and K, attempt to assign a reasonable value for theta for atom I
int AssignTheta(double&, int, int, int, Topology const&, Frame const&, std::vector<bool> const&, int);

} // END namespace Model
} // END namespace Structure
} // END namespace Cpptraj
#endif
