#ifndef INC_PUBFFT_H
#define INC_PUBFFT_H
/// C++ interface to Amber pubfft routines (pub_fft.F90)
class PubFFT {
  public:
    PubFFT();
    ~PubFFT();
    PubFFT(int);
    PubFFT(const PubFFT&);
    PubFFT& operator=(const PubFFT&);

    void Forward(double*);
    void Back(double*);
    int SetupFFT(int);
    int SetupFFTforN(int);
    void CorF_FFT(int, double*, double*);

    int size() { return fft_size_; }
  private:
    int fft_size_; ///< dimension of the FFT
    int saved_factors_size_;
    int saved_work_size_;
    int* saved_factors_;
    double* saved_work_;

    void Allocate();
};
#endif
