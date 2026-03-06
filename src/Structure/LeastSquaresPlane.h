#ifndef INC_STRUCTURE_LEASTSQUARESPLANE_H
#define INC_STRUCTURE_LEASTSQUARESPLANE_H
#include <vector>
#include "../Vec3.h"
class Frame;
class AtomMask;
class Vec3;
namespace Cpptraj {
namespace Structure {
/// Used to calculate the vector normal to a plane passing through a series of points
class LeastSquaresPlane {
  public:
    LeastSquaresPlane();
    /// Calculate vector normal to plane passing through atoms in given mask
    void CalcLeastSquaresPlane(Frame const&, AtomMask const&, bool);
    /// \return Vector normal to plane
    Vec3 const& Nxyz() const { return nxyz_; }
    /// \return Origin of vector normal to plane (plane center)
    Vec3 const& Cxyz() const { return cxyz_; }
  private:
    static double solve_cubic_eq(double,double,double,double);
    static Vec3 leastSquaresPlane(unsigned int, const double*);

    std::vector<double> vcorr_; ///< Will hold selected coordinates after centering.
    Vec3 nxyz_;                 ///< Hold vector normal to the plane
    Vec3 cxyz_;                 ///< Hold center (origin of vector normal to the plane)
};
}
}
#endif
