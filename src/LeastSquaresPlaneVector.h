#ifndef INC_LEASTSQUARESPLANEVECTOR_H
#define INC_LEASTSQUARESPLANEVECTOR_H
#include <vector>
#include "Vec3.h"
class Frame;
class AtomMask;
/// Calculate vector perpendicular to the plane passing through selected atoms.
class LeastSquaresPlaneVector {
  public:
    LeastSquaresPlaneVector();
    /// Reserve space for specified number of selected atoms.
    void ReserveForNumAtoms(unsigned int);
    /// \return vector perpendicular to the plane passing through selected atoms.
    Vec3 CalcLSPvec(Frame const&, AtomMask const&);
  private:
    static inline double solve_cubic_eq(double, double, double, double);
    static inline Vec3 leastSquaresPlane(int, const double*);

    std::vector<double> vcorr_;
};
#endif
