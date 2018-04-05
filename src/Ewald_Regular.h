#ifndef INC_EWALD_REGULAR_H
#define INC_EWALD_REGULAR_H
#include "Ewald.h"
/// Calculate regular Ewald energy
class Ewald_Regular : public Ewald {
  public:
    Ewald_Regular();
    /// Box, cutoff, dsum tol, rsum tol, ew coeff, maxexp, nb skin, erfc dx, debug, mlimits
    int Init(Box const&, double, double, double, double, double, double,
             double, int, const int*);
    // ----- Inherited ---------------------------
    int Setup(Topology const&, AtomMask const&);
    double CalcEnergy(Frame const&, AtomMask const&, double&); // TODO const?
  private:
    /// Determine max length for reciprocal calcs based on reciprocal limits
    static double FindMaxexpFromMlim(const int*, Matrix_3x3 const&);
    /// Determine max length for reciprocal calcs based on Ewald coefficient and recip tol.
    static double FindMaxexpFromTol(double, double);
    /// Determine reciprocal limits based on unit cell reciprocal vectors
    static void GetMlimits(int*, double, double, Vec3 const&, Matrix_3x3 const&);
    /// Ewald reciprocal energy
    double Recip_Regular(Matrix_3x3 const&, double);

    typedef Ewald::Darray Darray;
    // Hold trig tables
    Darray cosf1_;
    Darray cosf2_;
    Darray cosf3_;
    Darray sinf1_;
    Darray sinf2_;
    Darray sinf3_;
    Darray c12_;
    Darray s12_;
    Darray c3_;
    Darray s3_;
#   ifdef _OPENMP
    typedef Ewald::Iarray Iarray;
    Iarray mlim1_;        ///< Hold m1 reciprocal indices
    Iarray mlim2_;        ///< Hold m2 reciprocal indices
    int multCut_;         ///< Hold index after which multiplier should be 2.0.
#   endif
    double maxexp_;       ///< Determines how far out recip vectors go? FIXME check!
    double rsumTol_;      ///< Reciprocal space sum tolerance.
    int mlimit_[3];       ///< Number of units in each direction to calc recip. sum. / nfft
    int maxmlim_;         ///< The max of the three mlimit_ values. / pme spline order
};    
#endif
