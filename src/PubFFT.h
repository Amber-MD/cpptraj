#ifndef INC_PUBFFT_H
#define INC_PUBFFT_H
/// C++ interface to Amber pubfft routines (pubfft.F90)
class PubFFT {
  public:
    PubFFT();
    ~PubFFT();
    PubFFT(int);
    PubFFT(const PubFFT&);
    PubFFT& operator=(const PubFFT&);

    void Forward(double*);
    void Back(double*);
  private:
    int fft_size_; ///< dimension of the FFT
    int saved_factors_size_;
    int saved_work_size_;
    int* saved_factors_;
    double* saved_work_;
};
#endif
