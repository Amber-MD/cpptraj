#include <algorithm>
#include <cmath> // log, ldexp
#include "PubFFT.h"
#include "CpptrajStdio.h"

#ifndef FFTW_FFT
extern "C" {
  // pub_fft.F90 functions
  void pubfft_init_(int&, double*, int*);          // FFT init
  void pubfft_forward_(int&, double*, double*, int*); // Forward FFT
  void pubfft_back_(int&, double*, double*, int*); // Backward FFT
}
#endif

// CONSTRUCTOR
PubFFT::PubFFT() :
  fft_size_(0),
# ifdef FFTW_FFT
  in_(0)
# else
  saved_work_size_(0),
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
  if (saved_work_ != 0) delete[] saved_work_;
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
  saved_work_size_(rhs.saved_work_size_),
  saved_work_(0)
{
  std::copy(rhs.saved_factors_, rhs.saved_factors_+saved_factors_size_, saved_factors_);
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    std::copy(rhs.saved_work_, rhs.saved_work_+saved_work_size_, saved_work_);
  }
}
# endif

// Assignment operator
PubFFT& PubFFT::operator=(const PubFFT& rhs) {
  if (this == &rhs) return *this;
# ifdef FFTW_FFT
  // Nothing to copy, just allocate
  Allocate( rhs.fft_size_);
# else
  fft_size_ = rhs.fft_size_;
  if (saved_work_ != 0) delete[] saved_work_;
  std::copy(rhs.saved_factors_, rhs.saved_factors_+saved_factors_size_, saved_factors_);
  saved_work_size_ = rhs.saved_work_size_;
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    std::copy(rhs.saved_work_, rhs.saved_work_+saved_work_size_, saved_work_);
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
  pubfft_forward_( fft_size_, fft_array.CAptr(), saved_work_, saved_factors_ );
# endif
}

// PubFFT::Back()
void PubFFT::Back(ComplexArray& fft_array) {
# ifdef FFTW_FFT
  ComplexArrayToFFTW(fft_array);
  fftw_execute( backwards_plan_ );
  FFTWtoComplexArray(fft_array);
# else
  pubfft_back_( fft_size_, fft_array.CAptr(), saved_work_, saved_factors_);
# endif
}

// PubFFT::Allocate()
int PubFFT::Allocate(int sizeIn) {
  if (sizeIn < 0) {
    mprinterr("Error: Invalid memory size given for FFT (%i)\n", sizeIn);
    return 1;
  }
  fft_size_ = sizeIn;
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
  // NOTE: saved_factors_size is constant 
  if (saved_work_ != 0) delete[] saved_work_;
  // Reallocate FFT workspace
  std::fill(saved_factors_, saved_factors_+saved_factors_size_, 0);
  saved_work_size_ = 4 * fft_size_;
  if (saved_work_size_ > 0) {
    saved_work_ = new double[ saved_work_size_ ];
    std::fill(saved_work_, saved_work_+saved_work_size_, 0.0);
  } else if (saved_work_size_ < 0) {
    mprinterr("Error: Could not allocate memory for FFT; invalid size (%i)\n", saved_work_size_);
    return 1;
  } else
    saved_work_ = 0;
  // NOTE: Should this be called if fft_size is 0?
  pubfft_init_( fft_size_, saved_work_, saved_factors_ );
# endif
  return 0;
}

/** \return Next power of 2 >= sizeIn.
  * NOTE: This function is used when correlation is being calculated with
  *       FFT. As such the value is multiplied by 2.0 (via the final '+ 1'
  *       to provide enough space for zero padding to avoid end effects.
  */
inline static int NextPowOf2(int sizeIn) {
  return (ldexp( 1.0, (int)((log((double)sizeIn-1) / log(2.0)) + 1.0) + 1 ));
}

// PubFFT::SetupFFT_NextPowerOf2()
int PubFFT::SetupFFT_NextPowerOf2(int sizeIn) { return Allocate( NextPowOf2( sizeIn ) ); }

// PubFFT::SetupFFTforN()
int PubFFT::SetupFFTforN(int sizeIn) { return Allocate( sizeIn ); }
