#ifndef INC_SPLINE_H
#define INC_SPLINE_H
#include <vector>
/// Approximate curve using splines. 
class Spline {
  public:
    typedef std::vector<double> Darray;
    Spline() {}
    /// Calculate cubic spline coefficients for given X and Y values.
    int CubicSpline_Coeff(Darray const&, Darray const&);
    /// Calculate interpolated Y value from original X and Y values for given new X value.
    double CubicSpline_Eval(Darray const&, Darray const&, double) const;
    /// Calculate interpolated curve from original X and Y values for given new X values.
    Darray CubicSpline_Eval(Darray const&, Darray const&, Darray const&) const;
    /// Calculate interpolated Y value from regularly spaced original X values.
    inline double CubicSpline_Eval(double, double, Darray const&, double) const;
  private:
    // Cubic spline coefficients.
    Darray b_;
    Darray c_;
    Darray d_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
/** This version of Eval can be used when the original input X values
  * are regularly spaced.
  */
double Spline::CubicSpline_Eval(double del, double one_over_del, Darray const& yIn, double U) const
{
  int xidx = (int)(one_over_del * U);
  double dx = U - ((double)xidx * del);
  return yIn[xidx] + dx*(b_[xidx] + dx*(c_[xidx] + dx*d_[xidx]));
} 
#endif
