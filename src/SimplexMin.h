#ifndef INC_SIMPLEXMIN_H
#define INC_SIMPLEXMIN_H
#include "Random.h"
#include "DataSet.h"
class SimplexMin {
  public:
    typedef std::vector<double> Darray;
    /** Function that takes X values and parameters and calculates new Y values. */
    typedef double (*SimplexFunctionType)(DataSet*, Darray const&, Darray&);
    SimplexMin() : fxn_(0), Xvals_(0) {}
    int Minimize(SimplexFunctionType, Darray&, DataSet*, Darray const&,
                 double, int, double, Random_Number&);
  private:
    typedef std::vector<double>::size_type dsize;

    double chi_squared(Darray const& Ysearch);
    double Amotry(Darray& psum, int ihi, double fac);
    int Amoeba(int, double);
    void Average_vertices(Darray& xsearch) const;
    
    dsize NP_;                ///< Number of dimensions
    dsize NP1_;               ///< Number of vertices
    dsize Nvals_;             ///< Number of values being fit
    SimplexFunctionType fxn_; ///< Function to fit
    DataSet* Xvals_;          ///< DataSet containing X values (Nvals).
    Darray Xsimplex_;         ///< NP+1 rows by NP cols array containing current simplex.
    Darray Yvals_;            ///< Original y values (Nvals)
    Darray Ynew_;             ///< Y values with current parameters (Nvals)
    Darray Ysearch_;          ///< (NP1)
};
#endif
