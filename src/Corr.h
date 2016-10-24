#ifndef INC_CORR_H
#define INC_CORR_H
#include <vector>
#include "PubFFT.h"
/*! \file Corr.h
    \brief Classes that can be used to calculate time correlations from complex arrays.
 */
/// Used to directly calculate auto/cross-correlation for complex arrays.
/** Calculates correlation functions using the "direct" approach
  * (s. Comp. Sim. of Liquids, p.185)
  * - the result is not yet normalized by (no_of_discrete_data - t)**-1 (!)
  */
class CorrF_Direct {
  public:
    CorrF_Direct() : nsteps_(0) {}
    int CorrSetup(int);
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
    int CorrSetup(int);
    void AutoCorr(ComplexArray&);
    void CrossCorr(ComplexArray&, ComplexArray&);
    ComplexArray Array() { return ComplexArray( pubfft_.size() ); }
  private:
    PubFFT pubfft_;
};
#endif
