#include <cstdlib> // abs, intel 11 compilers choke on std::abs
#include <cmath> // sqrt
#include "DataSet_Vector.h"
#include "Constants.h" // For spherical harmonics norm.
#include "Corr.h"

const Vec3 DataSet_Vector::ZERO = Vec3(0,0,0);
const ComplexArray DataSet_Vector::COMPLEXBLANK = ComplexArray(0);

// CONSTRUCTOR
DataSet_Vector::DataSet_Vector() :
  DataSet(VECTOR, GENERIC, TextFormat(TextFormat::DOUBLE, 8, 4, 6), 1),
  order_(0) {}

size_t DataSet_Vector::MemUsageInBytes() const {
  size_t mySize = (vectors_.size() * Vec3::DataSize()) +
                  (origins_.size() * Vec3::DataSize()) +
                  sizeof(int);
  for (std::vector<ComplexArray>::const_iterator SH = sphericalHarmonics_.begin();
                                                 SH != sphericalHarmonics_.end();
                                               ++SH)
    mySize += (SH->DataSize());
  return mySize;
}
         

// DataSet_Vector::Allocate()
int DataSet_Vector::Allocate(SizeArray const& Nin) {
  if (!Nin.empty()) {
    vectors_.reserve( Nin[0] );
    origins_.reserve( Nin[0] ); // TODO: check if this needs allocation
  }
  return 0;
}

// DataSet_Vector::WriteBuffer()
void DataSet_Vector::WriteBuffer(CpptrajFile &cbuffer, SizeArray const& pIn) const {
  if (pIn[0] >= vectors_.size()) {
    cbuffer.Printf(format_.fmt(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); // VXYZ OXYZ
  } else {
    Vec3 const& Vxyz = vectors_[pIn[0]];
    Vec3 const& Oxyz = OXYZ(pIn[0]);
    cbuffer.Printf(format_.fmt(), Vxyz[0], Vxyz[1], Vxyz[2],
                                  Oxyz[0], Oxyz[1], Oxyz[2]);
  }
}

int DataSet_Vector::Append(DataSet* dsIn) {
  if (dsIn->Empty()) return 0;
  if (dsIn->Type() != VECTOR) return 1;
  Varray const& vIn = ((DataSet_Vector*)dsIn)->vectors_;
  Varray const& oIn = ((DataSet_Vector*)dsIn)->origins_;
  size_t oldsize = vectors_.size();
  vectors_.resize( oldsize + vIn.size() );
  std::copy( vIn.begin(), vIn.end(), vectors_.begin() + oldsize );
  if (oIn.empty() && !origins_.empty()) // Need vIn.size empty origin vectors
    origins_.resize( oldsize + vIn.size(), ZERO );
  else if (!oIn.empty() && origins_.empty()) // Need vectors_.size empty origin vecs
    origins_.resize( vectors_.size(), ZERO );
  if (!oIn.empty()) {
    oldsize = origins_.size();
    origins_.resize( oldsize + oIn.size() );
    std::copy( oIn.begin(), oIn.end(), origins_.begin() + oldsize );
  }
  return 0;
}

// DataSet_Vector::reset()
void DataSet_Vector::reset() {
  vectors_.clear();
  origins_.clear();
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

#ifdef MPI
int DataSet_Vector::Sync(size_t total, std::vector<int> const& rank_frames,
                         Parallel::Comm const& commIn)
{
  if (commIn.Size()==1) return 0;
  double buf[6];
  // TODO: Consolidate to 1 send/recv via arrays?
  if (commIn.Master()) {
    // Resize to accept data from other ranks.
    vectors_.resize( total );
    if (!origins_.empty()) {
      origins_.resize( total );
      int vidx = rank_frames[0]; // Index on master
      for (int rank = 1; rank < commIn.Size(); rank++) {
        for (int ridx = 0; ridx != rank_frames[rank]; ridx++, vidx++) {
          commIn.SendMaster( buf, 6, rank, MPI_DOUBLE );
          std::copy( buf,   buf+3, vectors_[vidx].Dptr() );
          std::copy( buf+3, buf+6, origins_[vidx].Dptr() );
        }
      }
    } else {
      int vidx = rank_frames[0]; // Index on master
      for (int rank = 1; rank < commIn.Size(); rank++) {
        for (int ridx = 0; ridx != rank_frames[rank]; ridx++, vidx++)
          commIn.SendMaster( vectors_[vidx].Dptr(), 3, rank, MPI_DOUBLE );
      } 
    }
  } else { // Send data to master
    if (!origins_.empty()) {
      for (unsigned int ridx = 0; ridx != vectors_.size(); ++ridx) {
        std::copy( vectors_[ridx].Dptr(), vectors_[ridx].Dptr()+3, buf   );
        std::copy( origins_[ridx].Dptr(), origins_[ridx].Dptr()+3, buf+3 );
        commIn.SendMaster( buf, 6, commIn.Rank(), MPI_DOUBLE );
      }
    } else
      for (unsigned int ridx = 0; ridx != vectors_.size(); ++ridx)
        commIn.SendMaster( vectors_[ridx].Dptr(), 3, commIn.Rank(), MPI_DOUBLE );
  }
  return 0;
}
#endif
