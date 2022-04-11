#ifndef INC_SPLINE_H
#define INC_SPLINE_H
#include <vector>
#include <cstddef> // size_t
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
    /// \return B coefficient array.
    Darray const& B_coeff() const { return b_; }
    /// \return C coefficient array.
    Darray const& C_coeff() const { return c_; }
    /// \return D coefficient array.
    Darray const& D_coeff() const { return d_; }
    /// /return memory usage in bytes
    size_t DataSize() const { return (b_.size() + c_.size() + d_.size()) * sizeof(double); }
  private:
    // Cubic spline coefficients.
    Darray b_;
    Darray c_;
    Darray d_;
};
#endif
