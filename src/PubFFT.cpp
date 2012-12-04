#include <cstring> // memset, memcpy
#include <cmath> // log, ldexp
#include "PubFFT.h"

extern "C" {
  // pub_fft.F90 functions
  void cffti_(int&, double*, int*);          // FFT init
  void cfftf_(int&, double*, double*, int*); // Forward FFT
  void cfftb_(int&, double*, double*, int*); // Backward FFT
}

inline static int NextPowOf2(int sizeIn) {
  return (ldexp( 1.0, (int)(log((double)(4 * sizeIn - 1)) / log(2.0)) + 1));
}

// CONSTRUCTOR
PubFFT::PubFFT() :
  fft_size_(0),
  saved_factors_size_(30),
  saved_work_size_(0),
  saved_factors_(0),
  saved_work_(0)
{}

// DESTRUCTOR
PubFFT::~PubFFT() {
  if (saved_factors_ != 0) delete[] saved_factors_;
  if (saved_work_ != 0) delete[] saved_work_;
}

// CONSTRUCTOR
/// Takes FFT size as input; ensures size is power of 2
PubFFT::PubFFT(int fft_sizeIn) :
  saved_factors_size_(30)
{
  int ndata = NextPowOf2( fft_sizeIn ); 
  fft_size_ = ndata / 2;
  saved_work_size_ = 2 * ndata; // 4 * fft_size
  saved_factors_ = new int[ saved_factors_size_ ];
  memset(saved_factors_, 0, saved_factors_size_ * sizeof(int));
  saved_work_ = new double[ saved_work_size_ ];
  memset(saved_work_, 0, saved_work_size_ * sizeof(double));
  cffti_( fft_size_, saved_work_, saved_factors_ );
}

// COPY CONSTRUCTOR
PubFFT::PubFFT(const PubFFT& rhs) :
  fft_size_(rhs.fft_size_),
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

// Assignment operator
PubFFT& PubFFT::operator=(const PubFFT& rhs) {
  if (this == &rhs) return *this;
  delete[] saved_work_;
  delete[] saved_factors_;
  fft_size_ = rhs.fft_size_;
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
  return *this;
}

// PubFFT::Forward()
void PubFFT::Forward(double* fft_array) {
  cfftf_( fft_size_, fft_array, saved_work_, saved_factors_ );
}

// PubFFT::Back()
void PubFFT::Back(double* fft_array) {
  cfftb_( fft_size_, fft_array, saved_work_, saved_factors_);
}

// PubFFT::Allocate()
void PubFFT::Allocate() {
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
}

// PubFFT::SetupFFT()
int PubFFT::SetupFFT(int sizeIn) {
  int ndata = NextPowOf2( sizeIn ); 
  fft_size_ = ndata / 2;
  saved_work_size_ = 2 * ndata;
  Allocate();
  return 0;
}

// PubFFT::SetupFFTforN()
int PubFFT::SetupFFTforN(int sizeIn) {
  fft_size_ = sizeIn;
  saved_work_size_ = 4 * fft_size_;
  Allocate();
  return 0;
}

void PubFFT::CorF_Auto(ComplexArray& data1) {
  cfftf_(fft_size_, data1.CAptr(), saved_work_, saved_factors_);
  // Calculate square modulus of F(data1)
  data1.SquareModulus();
  // Inverse FFT
  cfftb_(fft_size_, data1.CAptr(), saved_work_, saved_factors_);
  // Normalize with fft_size (since not done in inverse FFT routine)
  data1.Normalize( 1.0 / ((double) (fft_size_)));
}

void PubFFT::CorF_Cross(ComplexArray& data1, ComplexArray& data2) {
  // Cross-correlation
  cfftf_(fft_size_, data1.CAptr(), saved_work_, saved_factors_);
  cfftf_(fft_size_, data2.CAptr(), saved_work_, saved_factors_);
  // Calculate [data1]* x [data2] where * denotes complex conjugate.
  data1.ComplexConjTimes(data2);
  // Inverse FFT
  cfftb_(fft_size_, data1.CAptr(), saved_work_, saved_factors_);
  // Normalize with fft_size (since not done in inverse FFT routine)
  data1.Normalize( 1.0 / ((double) (fft_size_)));
}
