#ifndef INC_SPLINE_H
#define INC_SPLINE_H
#include <vector>
/// Approximate curve using splines. 
class Spline {
  public:
    typedef std::vector<double> Darray;
    Spline() {}
    void cubicSpline_coeff(Darray const&, Darray const&);
    Darray cubicSpline_eval(Darray const&, Darray const&, Darray const&) const;
  private:
    // Cubic spline coefficients.
    Darray b_;
    Darray c_;
    Darray d_;
};
#endif
