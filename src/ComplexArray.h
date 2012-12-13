#ifndef INC_COMPLEXARRAY_H
#define INC_COMPLEXARRAY_H
#include <vector>
#include "PubFFT.h"
/// Array that will hold complex numbers.
/** Implementation does not use STL vector so as to easily interface with
  * fortran routines for e.g. calculating FFT; current standard does not
  * directly allow access interal STL vector representation (i.e. double*).
  */
class ComplexArray {
  public:
    ComplexArray() : data_(0), ndata_(0), ncomplex_(0) {}
    ComplexArray(int);
    ComplexArray(ComplexArray const&);
    ComplexArray& operator=(ComplexArray const&);
    ~ComplexArray();
    void Allocate(int);
    void PadWithZero(int);
    /// Multiply all entries in data by input 
    void Normalize(double);
    void SquareModulus();
    /// Calculate [this]* x [rhs] where * denotes complex conjugate.
    void ComplexConjTimes(ComplexArray const&);
    double* CAptr()  { return data_;     }
    int size() const { return ncomplex_; }
    double& operator[](int idx)             { return data_[idx]; }
    double const& operator[](int idx) const { return data_[idx]; }
  private:
    double* data_;
    int ndata_;    ///< Actual size of data
    int ncomplex_; ///< # of complex numbers in array (ndata / 2)
};

/// Used to directly calculate auto/cross-correlation for complex arrays.
/** Calculates correlation functions using the "direct" approach
  * (s. Comp. Sim. of Liquids, p.185)
  * - the result is not yet normalized by (no_of_discrete_data - t)**-1 (!)
  */
class CorrF_Direct {
  public:
    CorrF_Direct() : nsteps_(0) {}
    CorrF_Direct(int stepsIn) : nsteps_(stepsIn), table_(2*nsteps_, 0.0) {}
    void Allocate(int);
    void AutoCorr(ComplexArray&);
    void CrossCorr(ComplexArray&, ComplexArray const&);
  private:
    int nsteps_;
    std::vector<double> table_;
};

/// Used to calculate auto/cross-correlation for complex arrays with FFTs.
class CorrF_FFT {
  public:
    CorrF_FFT() {}
    CorrF_FFT(int stepsIn) : pubfft_( stepsIn ) {}
    void Allocate(int);
    void AutoCorr(ComplexArray&);
    void CrossCorr(ComplexArray&, ComplexArray&);
    ComplexArray Array() { return ComplexArray( pubfft_.size() ); }
  private:
    PubFFT pubfft_;
};
#endif
