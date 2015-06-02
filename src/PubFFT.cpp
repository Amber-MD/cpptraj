#include <cstring> // memset, memcpy
#include <cmath> // log, ldexp
#include "PubFFT.h"

#ifndef FFTW_FFT
extern "C" {
  // pub_fft.F90 functions
  void cffti_(int&, double*, int*);          // FFT init
  void cfftf_(int&, double*, double*, int*); // Forward FFT
  void cfftb_(int&, double*, double*, int*); // Backward FFT
}
#endif

inline static int NextPowOf2(int sizeIn) {
  return (ldexp( 1.0, (int)(log((double)(4 * sizeIn - 1)) / log(2.0)) + 1));
}

// CONSTRUCTOR
PubFFT::PubFFT() :
  fft_size_(0),
# ifdef FFTW_FFT
  in_(0)
# else
  saved_factors_size_(30),
  saved_work_size_(0),
  saved_factors_(0),
  saved_work_(0)
# endif
{}

// DESTRUCTOR
PubFFT::~PubFFT() {
# ifdef FFTW_FFT
  if (in_ != 0) {
    fftw_free(in_);
    fftw_destroy_plan(forwards_plan_);
    fftw_destroy_plan(backwards_plan_);
  }
# else
  if (saved_factors_ != 0) delete[] saved_factors_;
  if (saved_work_ != 0) delete[] saved_work_;
# endif
}

// CONSTRUCTOR
/// Takes FFT size as input; ensures size is power of 2
PubFFT::PubFFT(int fft_sizeIn)
# ifndef FFTW_FFT
  : saved_factors_size_(30)
# endif
{
  int ndata = NextPowOf2( fft_sizeIn ); 
  fft_size_ = ndata / 2;
# ifdef FFTW_FFT
  in_  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_size_);
  forwards_plan_  = fftw_plan_dft_1d(fft_size_, in_, in_, FFTW_FORWARD,  FFTW_ESTIMATE);
  backwards_plan_ = fftw_plan_dft_1d(fft_size_, in_, in_, FFTW_BACKWARD, FFTW_ESTIMATE);
# else
  saved_work_size_ = 2 * ndata; // 4 * fft_size
  saved_factors_ = new int[ saved_factors_size_ ];
  memset(saved_factors_, 0, saved_factors_size_ * sizeof(int));
  saved_work_ = new double[ saved_work_size_ ];
  memset(saved_work_, 0, saved_work_size_ * sizeof(double));
  cffti_( fft_size_, saved_work_, saved_factors_ );
# endif
}

// COPY CONSTRUCTOR
PubFFT::PubFFT(const PubFFT& rhs) :
  fft_size_(rhs.fft_size_),
# ifdef FFTW_FFT
  in_(0)
{
  if (fft_size_ > 0) {
    in_  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_size_);
    forwards_plan_  = fftw_plan_dft_1d(fft_size_, in_, in_, FFTW_FORWARD,  FFTW_ESTIMATE);
    backwards_plan_ = fftw_plan_dft_1d(fft_size_, in_, in_, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
}
# else
  saved_factors_size_(rhs.saved_factors_size_),
  saved_work_size_(rhs.saved_work_size_),
  saved_factors_(0),
  saved_work_(0)
{
  if (saved_factors_size_ > 0) {
    saved_factors_ = new int[ saved_factors_size_ ];
    memcpy(saved_factors_, rhs.saved_factors_, saved_factors_size_ * sizeof(int));
  }
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    memcpy(saved_work_, rhs.saved_work_, saved_work_size_ * sizeof(double));
  }
}
# endif

// Assignment operator
PubFFT& PubFFT::operator=(const PubFFT& rhs) {
  if (this == &rhs) return *this;
  fft_size_ = rhs.fft_size_;
# ifdef FFTW_FFT
  Allocate();
# else
  if (saved_work_ != 0) delete[] saved_work_;
  if (saved_factors_ != 0) delete[] saved_factors_;
  saved_factors_size_ = rhs.saved_factors_size_;
  saved_work_size_ = rhs.saved_work_size_;
  if (saved_factors_size_ > 0) {
    saved_factors_ = new int[ saved_factors_size_ ];
    memcpy(saved_factors_, rhs.saved_factors_, saved_factors_size_ * sizeof(int));
  } else
    saved_factors_ = 0;
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    memcpy(saved_work_, rhs.saved_work_, saved_work_size_ * sizeof(double));
  } else
    saved_work_ = 0;
# endif
  return *this;
}

#ifdef FFTW_FFT
void PubFFT::ComplexArrayToFFTW(ComplexArray const& fft_array) {
  // TODO: Can we just use the ComplexArray pointer or direct copy?
  unsigned int idx0 = 0; // Index into fft_array
  for (int idx1 = 0; idx1 != fft_array.size(); idx1++, idx0 += 2) {
    in_[idx1][0] = fft_array[idx0  ];
    in_[idx1][1] = fft_array[idx0+1];
  }
}

void PubFFT::FFTWtoComplexArray(ComplexArray& fft_array) const {
  unsigned int idx0 = 0; // Index into fft_array
  for (int idx1 = 0; idx1 != fft_array.size(); idx1++, idx0 += 2) {
    fft_array[idx0  ] = in_[idx1][0];
    fft_array[idx0+1] = in_[idx1][1];
  }
}
#endif

// PubFFT::Forward()
void PubFFT::Forward(ComplexArray& fft_array) {
# ifdef FFTW_FFT
  ComplexArrayToFFTW(fft_array);
  fftw_execute( forwards_plan_ );
  FFTWtoComplexArray(fft_array);
# else
  cfftf_( fft_size_, fft_array.CAptr(), saved_work_, saved_factors_ );
# endif
}

// PubFFT::Back()
void PubFFT::Back(ComplexArray& fft_array) {
# ifdef FFTW_FFT
  ComplexArrayToFFTW(fft_array);
  fftw_execute( backwards_plan_ );
  FFTWtoComplexArray(fft_array);
# else
  cfftb_( fft_size_, fft_array.CAptr(), saved_work_, saved_factors_);
# endif
}

// PubFFT::Allocate()
void PubFFT::Allocate() {
# ifdef FFTW_FFT
  // Destroy any existing plans/storage space, allocate new space.
  if (in_ != 0) {
    fftw_free(in_);
    fftw_destroy_plan(forwards_plan_);
    fftw_destroy_plan(backwards_plan_);
  }
  in_  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_size_);
  forwards_plan_  = fftw_plan_dft_1d(fft_size_, in_, in_, FFTW_FORWARD,  FFTW_ESTIMATE);
  backwards_plan_ = fftw_plan_dft_1d(fft_size_, in_, in_, FFTW_BACKWARD, FFTW_ESTIMATE);
# else
  // NOTE: Do not change saved_factors_size
  if (saved_work_ != 0) delete[] saved_work_;
  if (saved_factors_ != 0) delete[] saved_factors_;
  // Reallocate FFT workspace
  if (saved_factors_size_ > 0) {
    saved_factors_ = new int[ saved_factors_size_ ];
    memset(saved_factors_, 0, saved_factors_size_ * sizeof(int));
  } else
    saved_factors_ = 0;
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    memset(saved_work_, 0, saved_work_size_ * sizeof(double));
  } else
    saved_work_ = 0;
  // NOTE: Should this be called if fft_size is 0?
  cffti_( fft_size_, saved_work_, saved_factors_ );
# endif
}

// PubFFT::SetupFFT_NextPowerOf2()
int PubFFT::SetupFFT_NextPowerOf2(int sizeIn) {
  int ndata = NextPowOf2( sizeIn ); 
  fft_size_ = ndata / 2;
# ifndef FFTW_FFT
  saved_work_size_ = 2 * ndata;
# endif
  Allocate();
  return 0;
}

// PubFFT::SetupFFTforN()
int PubFFT::SetupFFTforN(int sizeIn) {
  fft_size_ = sizeIn;
# ifndef FFTW_FFT
  saved_work_size_ = 4 * fft_size_;
# endif
  Allocate();
  return 0;
}
