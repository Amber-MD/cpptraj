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
    PubFFT(int);
    PubFFT(const PubFFT&);
    PubFFT& operator=(const PubFFT&);

    int size() const { return fft_size_; }
    void Forward(ComplexArray&);
    void Back(ComplexArray&);
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
    int saved_factors_size_;
    int saved_work_size_;
    int* saved_factors_;
    double* saved_work_;
#   endif
    void Allocate();
};
#endif
