#ifndef INC_SIMPLEXMIN_H
#define INC_SIMPLEXMIN_H
#include <vector>
#include "Random.h"
class SimplexMin {
  public:
    typedef std::vector<double> Darray;
    /** Function that takes parameters can calculates new Y values. */
    typedef double (*SimplexFunctionType)(Darray const&, Darray&);
    SimplexMin();
  private:
    typedef std::vector<double>::size_type dsize;

    double chi_squared(Darray const& Ysearch);
    double Amotry(Darray& psum, int ihi, double fac);
    int Amoeba(int, double);
    void Average_vertices(Darray& xsearch) const;
    int Simplex_min(SimplexFunctionType fxnIn, Darray& Q_vector,
                            Darray const& YvalsIn, double delqfracIn, int, double,
                            Random_Number&);
    
    Darray Xsimplex_; ///< NP+1 rows by NP cols array containing current simplex.
    dsize NP_; ///< Number of dimensions
    dsize NP1_; ///< Number of vertices
    dsize Nvals_; ///< Number of values being fit
    SimplexFunctionType fxn_; ///< Function to fit
    Darray Yvals_; ///< Original y values (Nvals)
    Darray Ynew_; ///< Y values with current parameters (Nvals)
    Darray Ysearch_; // NP1
};
#endif
