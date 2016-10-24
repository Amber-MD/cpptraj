#ifndef INC_PUBFFT_H
#define INC_PUBFFT_H
#include "ComplexArray.h"
#ifdef FFTW_FFT
# include <fftw3.h>
#endif
/// C++ interface to Amber pubfft routines (pub_fft.F90)
class PubFFT {
  public:
    PubFFT();
    ~PubFFT();
    PubFFT(const PubFFT&);
    PubFFT& operator=(const PubFFT&);
    /// \return FFT size in terms of number of complex numbers.
    int size() const { return fft_size_; }
    void Forward(ComplexArray&);
    void Back(ComplexArray&);
    /// Set up FFT with size == to next power of 2, times 2 for zero padding
    int SetupFFT_NextPowerOf2(int);
    int SetupFFTforN(int);
  private:
    int fft_size_; ///< dimension of the FFT
#   ifdef FFTW_FFT
    void ComplexArrayToFFTW(ComplexArray const&);
    void FFTWtoComplexArray(ComplexArray&) const;

    fftw_complex* in_;
    fftw_plan forwards_plan_;
    fftw_plan backwards_plan_;
#   else
    static const int saved_factors_size_ = 30;
    int saved_work_size_;
    int saved_factors_[saved_factors_size_];
    double* saved_work_;
#   endif
    int Allocate(int);
};
#endif
