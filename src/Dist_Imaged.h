#ifndef INC_DIST_IMAGED_H
#define INC_DIST_IMAGED_H
class Vec3;
class Matrix_3x3;
namespace Cpptraj {
/// Calculate minimum imaged distance (squared) between two coordinates in Cartesian space
double Dist2_Imaged(Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&, int*);
}
#endif
