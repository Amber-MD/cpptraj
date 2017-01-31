#ifndef INC_SPLINE_H
#define INC_SPLINE_H
#include <vector>
/// Approximate curve using splines. 
class Spline {
  public:
    typedef std::vector<double> Darray;
    Spline() {}
    int CubicSpline_Coeff(Darray const&, Darray const&);
    Darray CubicSpline_Eval(Darray const&, Darray const&, Darray const&) const;
  private:
    // Cubic spline coefficients.
    Darray b_;
    Darray c_;
    Darray d_;
};
#endif
