#ifndef INC_STRUCTURE_MODEL_H
#define INC_STRUCTURE_MODEL_H
class Topology;
class Frame;
namespace Cpptraj {
namespace Structure {
/// Routines to generate model parameters
namespace Model {
/// Given atoms J K and L, attempt to assign a reasonable value for phi for atom I
int AssignPhi(double&, int, int, int, int, Topology const&, Frame const&);

} // END namespace Model
} // END namespace Structure
} // END namespace Cpptraj
#endif
