#ifndef INC_SPLINE_H
#define INC_SPLINE_H
#include <vector>
class Spline {
  public:
    typedef std::vector<double> Darray;
    Spline() {}
    /// Calculate cubic spline coefficients for given X and Y values.
    int CalcCoeff_Cubic(Darray const&, Darray const&);
    /// Given an X value and original X and Y values, output splined Y value.
    double SplinedY(double, Darray const&, Darray const&) const;
    /// Given an array of new X values and original X and Y values, output splined Y values.
    Darray SplinedYvals(Darray const&, Darray const&, Darray const&) const;
  private:
    Darray B_;
    Darray C_;
    Darray D_;
};
#endif
