#ifndef INC_CURVEFIT_H
#define INC_CURVEFIT_H
#ifdef DBG_CURVEFIT
# include <cstdio> // DEBUG
#endif
#include <vector>
/// Used for non-linear curve fitting.
class CurveFit {
  public:
    typedef std::vector<double> Darray;
    /** Function prototype for Y = f(X):
      * X value, P[] parameters, dYdP[] first derivatives of Y w.r.t. P[]
      * \return Y
      */
    //typedef double (*FitFunctionType)(double, Darray const&, Darray&);
    /** Function prototype for Y =f(x):
      * Xvalue, P[] parameters.
      */
    typedef double (*FitFunctionType)(double, Darray const&);
    
    CurveFit();
#   ifdef DBG_CURVEFIT
    ~CurveFit(); // DEBUG
#   endif
    /// Perform Levenberg-Marquardt curve fit: fxn, x, y, p, tol, iter
    int LevenbergMarquardt(FitFunctionType, Darray const&, Darray const&, Darray&,
                           double, int);
    /// \return Status message.
    static const char* Message(int);
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<double>::size_type dsize;
    static const double machine_epsilon;

    /// Calculate residual only
    void EvaluateFxn(Darray const&, Darray const&, Darray const&, Darray&) const;
    /// Calculate Jacobian using forward-difference approximation
    void CalcJacobian_ForwardDiff(Darray const&, Darray const&, Darray&, Darray const&, Darray&);
    /// Calculate || m(i,...) || for vector in matrix
    static double VecNorm(Darray::const_iterator const&, dsize);
    /// Calculate || v || for vector
    static inline double VecNorm(Darray const& vec) {
      return VecNorm( vec.begin(), vec.size() );
    }
    /// Print final parameters to STDOUT
    void PrintFinalParams(Darray const&) const;
    /// \return true if input coords/parameters have problems.
    bool ParametersHaveProblems(Darray const&, Darray const&, Darray const&) const;
    // DEBUG
    inline void PrintMatrix(const char*, int, int, Darray const&) const;
    inline void PrintVector(const char*, Darray const&) const;
    void DBGPRINT(const char*, ...) const;

    FitFunctionType fxn_; ///< Function to fit to.
    dsize m_; ///< Number of values (rows)
    dsize n_; ///< Number of parameters (cols)
    Darray jacobian_; ///< Jacobian/R, stored in transpose (row-major)
    // DEBUG
#   ifdef DBG_CURVEFIT
    FILE* dbgfile_;
#   endif
};
#endif
