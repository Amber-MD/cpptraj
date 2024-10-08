#ifndef INC_ENERGY_EWALD_RECIP_H
#define INC_ENERGY_EWALD_RECIP_H
#include <vector>
class Box;
class EwaldOptions;
class Matrix_3x3;
class Vec3;
namespace Cpptraj {
namespace Energy {
class EwaldParams;
/// For calculating reciprocal space energy using the "regular" Ewald summation
class Ewald_Recip {
  public:
    Ewald_Recip();


    int InitRecip(EwaldOptions const&, EwaldParams const&, Box const&, int);
    void PrintRecipOpts() const;
    /// Set up trig tables for the given number of selected atoms
    int SetupRecip(int);
  private:
    typedef std::vector<double> Darray;

    /// Determine max length for reciprocal calcs based on reciprocal limits
    static double FindMaxexpFromMlim(const int*, Matrix_3x3 const&);
    /// Determine max length for reciprocal calcs based on Ewald coefficient and recip tol.
    static double FindMaxexpFromTol(double, double);
    /// Determine reciprocal limits based on unit cell reciprocal vectors
    static void GetMlimits(int*, double, double, Vec3 const&, Matrix_3x3 const&);

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
    double maxexp_;       ///< Determines how far out recip vectors go? TODO check!
    double rsumTol_;      ///< Reciprocal space sum tolerance.
    int mlimit_[3];       ///< Number of units in each direction to calc recip. sum.
    int maxmlim_;         ///< The max of the three mlimit_ values.
    int debug_;
};
}
}
#endif
