#ifndef INC_SPLINEFXNTABLE_H
#define INC_SPLINEFXNTABLE_H
#include <vector>
//#incl ude "CpptrajStdio.h" // DEBUG
/// Can be used to approximate a function using cubic splines.
class SplineFxnTable {
  public:
    /// CONSTRUCTOR
    SplineFxnTable();
    /// Generic form of the function to approximate.
    typedef double (*FxnType)(double);
    /// Fill the ftable using given function and spacing, from given min to max.
    int FillTable(FxnType, double, double, double);
    /// Fill the table using given function, mesh size, from given min to max.
    int FillTable(FxnType, int, double, double);
    /// \return Approximated Y value from given X value.
    double Yval(double xIn) const {
      double xval = xIn - Xmin_;
      long int xidx = ((long int)(one_over_Dx_ * xval));
      //mprintf("SPLINEFXNTBL: idx= %li\n", xidx); // DEBUG
#     ifdef SPLINEFXNTABLE_CHECK_RANGE
      // Protect against out of range values
      // NOTE - If args to Yval are chosen carefully, omitting this check speeds things up.
      if (xidx >= (long int)Xvals_.size())
        xidx = Xvals_.size() - 1;
      else if (xidx < 0)
        xidx = 0;
#     endif
      // Delta from index
      double dx = xval - ((double)xidx * Dx_);
      // DEBUG
      //mprintf("DEBUG: xidx= %8li  x=%20.10E  x*dx= %20.10E  xval= %20.10E  dx= %20.10E  newDx= %20.10E\n",
      //        xidx, xIn, (double)(xidx)*Dx_, Xvals_[xidx], dx, xIn - Xvals_[xidx]);
      // Index into the table
      xidx *= 4;
      // DEBUG
      //if (xidx < 0 || xidx >= (int)table_.size())
      //  mprinterr("Error: index %li out of range (%zu) for X val %g (Xmin= %g 1/dx= %g)\n", xidx, table_.size(), xIn, Xmin_, one_over_Dx_);
      return table_[xidx] + 
             dx*(table_[xidx+1] + dx*(table_[xidx+2] + dx*table_[xidx+3]));
    }
    /// \return Approximated Y value, use internal X table
    double Yval_xtable(double xIn) const {
      long int xidx = (long int)(one_over_Dx_ * (xIn - Xmin_));
#     ifdef SPLINEFXNTABLE_CHECK_RANGE
      // Protect against out of range values
      // NOTE - If args to Yval are chosen carefully, omitting this check speeds things up.
      if (xidx >= (long int)Xvals_.size())
        xidx = Xvals_.size() - 1;
      else if (xidx < 0)
        xidx = 0;
#     endif
      // Delta from index
      double dx = xIn - Xvals_[xidx];
      // Index into the table
      xidx *= 4;
      return table_[xidx] + 
             dx*(table_[xidx+1] + dx*(table_[xidx+2] + dx*table_[xidx+3]));
    }
    /// \return More accurate X value via search of internal X table
    double Yval_accurate(double) const;
    /// \return Number of values in the table
    unsigned int Nvals() const { return Xvals_.size(); }
    /// Print table details to STDOUT
    void PrintTableInfo(const char*) const;
    /// Print memory usage to STDOUT
    void PrintMemUsage(const char*) const;
    // DEBUG Access to internal table
    std::vector<double> const& InternalTable() const { return table_; }
  private:
    typedef std::vector<double> Darray;

    double Dx_;          ///< Spacing
    double one_over_Dx_; ///< 1 over spacing
    double Xmin_;        ///< Minimum value for which function can be approximated
    double Xmax_;        ///< Maximum value for which function can be approximated
    Darray table_;       ///< Hold Y followed by spline B C D coefficients
    Darray Xvals_;       ///< Hold X values
};
#endif
