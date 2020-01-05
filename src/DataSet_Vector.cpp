#include <cstdlib> // abs, intel 11 compilers choke on std::abs
#include <cmath> // sqrt
#include "DataSet_Vector.h"
#include "Constants.h" // For spherical harmonics norm.
#include "Corr.h"

const ComplexArray DataSet_Vector::COMPLEXBLANK = ComplexArray(0);

/// CONSTRUCTOR
DataSet_Vector::DataSet_Vector() : order_(0) {}

// CONSTRUCTOR
DataSet_Vector::DataSet_Vector(DataType tIn, TextFormat const& fIn) :
  DataSet(tIn, VECTOR_1D, fIn, 1),
  order_(0) {}

/** \return size of vectors plus spherical harmonics. */
size_t DataSet_Vector::internalSize() const {
  size_t mySize = (vectors_.size() * Vec3::DataSize()) +
                  sizeof(int);
  for (std::vector<ComplexArray>::const_iterator SH = sphericalHarmonics_.begin();
                                                 SH != sphericalHarmonics_.end();
                                               ++SH)
    mySize += (SH->DataSize());
  return mySize;
}

/** Allocate vector array. */
void DataSet_Vector::internalAlloc(SizeArray const& Nin) {
  vectors_.reserve( Nin[0] );
}

/** Append to vector array. */
void DataSet_Vector::internalAppend(const DataSet_Vector* dsIn) {
  Varray const& vIn = dsIn->vectors();
  size_t oldsize = vectors_.size();
  vectors_.resize( oldsize + vIn.size() );
  std::copy( vIn.begin(), vIn.end(), vectors_.begin() + oldsize );
}

/** Reset the vector array */
void DataSet_Vector::internalReset() {
  vectors_.clear();
  sphericalHarmonics_.clear();
  order_ = 0;
}

// -----------------------------------------------------------------------------
int DataSet_Vector::CalcVectorCorr(DataSet_Vector const& V2, DataSet_1D& Ct,
                                   int lagmaxIn) const
{
  if (Ct.Type() != DataSet::DOUBLE)
    return 1;
  unsigned int Nvecs = this->Size();
  if (Nvecs != V2.Size()) return 1;
  if (Nvecs < 2) return 1;
  unsigned int lagmax;
  if (lagmaxIn == -1)
    lagmax = Nvecs;
  else if (lagmaxIn > (int)Nvecs)
    lagmax = Nvecs;
  else
    lagmax = (unsigned int)lagmaxIn;
  bool crosscorr = (&V2 != this);
  
  unsigned int arraySize = Nvecs * 3; // XYZ
  CorrF_FFT pubfft;
  pubfft.CorrSetup( arraySize );
  ComplexArray data1 = pubfft.Array();
  data1.PadWithZero( arraySize );
  ComplexArray data2;
  if (crosscorr)
    data2 = data1;
  
  // Load up real components of data1 with vector XYZ
  unsigned int idx = 0;
  for (unsigned int v = 0; v != Nvecs; v++, idx += 6) {
    data1[idx  ] = vectors_[v][0];
    data1[idx+2] = vectors_[v][1];
    data1[idx+4] = vectors_[v][2];
    if (crosscorr) {
      data2[idx  ] = V2.vectors_[v][0];
      data2[idx+2] = V2.vectors_[v][1];
      data2[idx+4] = V2.vectors_[v][2];
    }
  }

  if (crosscorr)
    pubfft.CrossCorr(data1, data2);
  else
    pubfft.AutoCorr(data1);

  // Place desired components of correlation fn back in output and normalize
  double normV = (double)arraySize;
  double norm0 = 1.0 / (fabs(data1[0]) / normV);
  idx = 0;
  for (unsigned int iout = 0; iout != lagmax; iout++, normV -= 3.0, idx += 6) {
    double c_t = (data1[idx] / normV) * norm0;
    Ct.Add(iout, &c_t);
    //Ct[iout] = (data1[idx] / normV) * norm0;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// DataSet_Vector::CalcSphericalHarmonics()
/** Calc spherical harmonics of order l=0,1,2 and -l<=m<=l for stored vectors 
  * (see e.g. Merzbacher, Quantum Mechanics, p. 186).
  */
int DataSet_Vector::CalcSphericalHarmonics(int orderIn) {
  const double SH00=0.28209479;
  const double SH10=0.48860251;
  const double SH11=0.34549415;
  const double SH20=0.31539157;
  const double SH21=0.77254840;
  const double SH22=0.38627420;
  if (orderIn < 0 || orderIn > 2) return 1;
  if (order_ == orderIn) return 0;
  order_ = orderIn;
  sphericalHarmonics_.clear();
  // Allocate and init complex arrays
  sphericalHarmonics_.resize( (order_ * 2) + 1, ComplexArray(Size()) );
  //Loop over all vectors, convert to spherical harmonics
  unsigned int cidx = 0; // Index into complex arrays
  for (Varray::const_iterator vec = vectors_.begin(); 
                              vec != vectors_.end(); ++vec, cidx += 2)
  {
    const Vec3& Vxyz = *vec;
    double len = sqrt( Vxyz.Magnitude2() );
    // Loop over m = -legendre order, ..., +legendre_order
    for (int midx = -order_; midx <= order_; ++midx) {
      ComplexArray& D = sphericalHarmonics_[midx + order_];
       // D[cidx] is the real part, D[cidx+1] is the imaginary part.
      double ri = 1.0 / len;
      if (order_ == 0 && midx == 0) {
        D[cidx] = SH00;
      } else if (order_ == 1) {
        if (midx == 0) {
          D[cidx] = SH10 * Vxyz[2] * ri;
        } else {
          D[cidx  ] = -midx * SH11 * Vxyz[0] * ri;
          D[cidx+1] =        -SH11 * Vxyz[1] * ri;
        }
      } else if (order_ == 2) {
        if (midx == 0) {
          D[cidx] = SH20 * (2.0*Vxyz[2]*Vxyz[2] - Vxyz[0]*Vxyz[0] - Vxyz[1]*Vxyz[1]) * ri * ri;
        } else if (abs(midx) == 1) {
          D[cidx  ] = -midx * SH21 * Vxyz[0] * Vxyz[2] * ri * ri;
          D[cidx+1] =        -SH21 * Vxyz[1] * Vxyz[2] * ri * ri;
        } else {
          D[cidx  ] = SH22 * (Vxyz[0]*Vxyz[0] - Vxyz[1]*Vxyz[1]) * ri * ri;
          D[cidx+1] = midx * SH22 * Vxyz[0] * Vxyz[1] * ri * ri;
        }
      }
      //mprintf("DBG: Vec %zu sphereHarm(m=%i) = %f + %fi\n",vec-vectors_.begin(), midx,
      //        D[cidx], D[cidx+1]);
    }
  }
  return 0; 
}

// DataSet_Vector::SphericalHarmonics()
ComplexArray const& DataSet_Vector::SphericalHarmonics(int midx) const {
  if (sphericalHarmonics_.empty() || abs(midx) > order_)
   return COMPLEXBLANK;
  return sphericalHarmonics_[midx + order_];
}

/** 4*PI / ((2*order)+1) due to spherical harmonics addition theorem */
double DataSet_Vector::SphericalHarmonicsNorm(int order) {
  if      (order == 2) return Constants::FOURFIFTHSPI;
  else if (order == 1) return Constants::FOURTHIRDSPI;
  else if (order == 0) return Constants::FOURPI;
  else return 1.0;
}
