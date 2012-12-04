#include <cstring> // memcpy
#include <cstdlib> // abs, intel 11 compilers choke on std::abs
#include <cmath> // sqrt
#include "DataSet_Vector.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_Vector::DataSet_Vector() :
  DataSet(VECTOR, 8, 4, 1),
  totalidx_(0),
  currentidx_(0),
  order_(0),
  xyz_(0),
  writeSum_(false),
  isIred_(false)
{ }

// DESTRUCTOR
DataSet_Vector::~DataSet_Vector() {
  if (xyz_!=0) delete[] xyz_;
}

// DataSet_Vector::IncreaseSize()
void DataSet_Vector::IncreaseSize() {
  int newsize = totalidx_ + 3000; // 500 frames * 6
  double* newxyz = new double[ newsize ];
  if (totalidx_ > 0) {
    memcpy(newxyz,        xyz_, totalidx_ * sizeof(double));
    delete[] xyz_;
  }
  memset(newxyz+totalidx_, 0,      3000 * sizeof(double));
  totalidx_ = newsize;
  xyz_ = newxyz;
}

// DataSet_Vector::Allocate()  
int DataSet_Vector::Allocate(int Nin) {
  if (xyz_!=0) delete[] xyz_;
  if (Nin < 1) {
    // # of frames is not known. Allocation will happen via IncreaseSize( )
    totalidx_ = 0;
    xyz_ = 0;
  } else {
    totalidx_ = 6 * Nin;
    xyz_ = new double[ totalidx_ ];
    memset( xyz_, 0, totalidx_ * sizeof(double) );
  }
  currentidx_ = 0;
    
  return 0;
}

// DataSet_Vector::WriteBuffer()
void DataSet_Vector::WriteBuffer(CpptrajFile &cbuffer, int frameIn) {
  int idx = frameIn * 6;
  if (idx < 0 || frameIn >= currentidx_) {
    mprinterr("Error: DataSet_Vector: Frame %i is out of range.\n",frameIn);
    return;
  }
  cbuffer.Printf(data_format_, xyz_[idx  ]);
  cbuffer.Printf(data_format_, xyz_[idx+1]);
  cbuffer.Printf(data_format_, xyz_[idx+2]);
  cbuffer.Printf(data_format_, xyz_[idx+3]);
  cbuffer.Printf(data_format_, xyz_[idx+4]);
  cbuffer.Printf(data_format_, xyz_[idx+5]);
  if (writeSum_) {
    cbuffer.Printf(data_format_, xyz_[idx  ] + xyz_[idx+3]);
    cbuffer.Printf(data_format_, xyz_[idx+1] + xyz_[idx+4]);
    cbuffer.Printf(data_format_, xyz_[idx+2] + xyz_[idx+5]);
  }
}

// DataSet_Vector::sphericalHarmonics()
/** Calc spherical harmonics of order l=0,1,2
  * and -l<=m<=l with cartesian coordinates as input
  * (see e.g. Merzbacher, Quantum Mechanics, p. 186)
  * D[0] is the real part, D[1] is the imaginary part.
  */
void DataSet_Vector::sphericalHarmonics(int l, int m, const double* XYZ, double r, double D[2])
{
  const double SH00=0.28209479;
  const double SH10=0.48860251;
  const double SH11=0.34549415;
  const double SH20=0.31539157;
  const double SH21=0.77254840;
  const double SH22=0.38627420;

  double ri;
  double x = XYZ[0];
  double y = XYZ[1];
  double z = XYZ[2];
  //mprintf("CDBG: x=%.10lE y=%.10lE z=%.10lE\n",x,y,z);

  D[0] = 0.0;
  D[1] = 0.0;
  ri = 1.0 / r;

  if (l == 0 && m == 0) {
    D[0] = SH00;
  } else if (l == 1) {
    if (m == 0) {
      D[0] = SH10 * z * ri;
    } else {
      D[0] = -m * SH11 * x * ri;
      D[1] =     -SH11 * y * ri;
    }
  } else if(l == 2) {
    if (m == 0) {
      D[0] = SH20 * (2.0*z*z - x*x - y*y) * ri * ri;
    } else if (abs(m) == 1) {
      D[0] = -m * SH21 * x * z * ri * ri;
      D[1] =     -SH21 * y * z * ri * ri;
    } else {
      D[0] = SH22 * (x*x - y*y) * ri * ri;
      D[1] = m * SH22 * x * y * ri * ri;
    }
  }
  //mprintf("CDBG: dreal=%.10lE dimg=%.10lE\n",D[0],D[1]);
}

// DataSet_Vector::CalcSphericalHarmonics()
/** Calculate spherical harmonics for all vectors in DataSet for the
  * specified order. Spherical harmonics will be stored internally
  * in the sphereharm array, which will speed up any subsequent calcs
  * at the expense of memory. The layout is currently e.g. (order 2):
  *   V0[m-2R], V0[m-2C], V0[m-1R], V0[m-1C], ... V0[m+2C], V1[m-2R], ...
  */
// TODO: Just return an array?
// TODO: Change memory layout so all values for a given order for each 
//       vector follow each other, i.e. order 2: 
//    V0[m-2R], V0[m-2C], ... VN[m-2R], VN[m-2C], V0[m-1R], V0[m-1C], ...
void DataSet_Vector::CalcSphericalHarmonics(int orderIn) {
  double Dcomplex[2];
  // If sphereharm is not empty assume SphericalHarmonics was already called.
  if (!sphereharm_.empty()) return;
  // Store order for checking calls to FillData
  order_ = orderIn;
  // Allocate space to hold complex numbers. Each frame has 2 doubles
  // (real + imaginary)  for each m = -olegendre -> +olegendre
  int n_of_vecs = Size();
  int p2blocksize = 2 * ((order_ * 2) + 1);
  int p2size = (p2blocksize * n_of_vecs);
  sphereharm_.reserve( p2size );
   // Loop over all vectors, convert to spherical harmonics
  double* VXYZ = xyz_;
  // TODO: Use iterator?
  for (int i = 0; i < n_of_vecs; ++i) {
    // Calc vector length
    double len = sqrt(VXYZ[0]*VXYZ[0] + VXYZ[1]*VXYZ[1] + VXYZ[2]*VXYZ[2]);
    // Loop over m = -olegendre, ..., +olegendre
    for (int midx = -order_; midx <= order_; ++midx) {
      sphericalHarmonics(order_, midx, VXYZ, len, Dcomplex);
      //mprinterr("DBG: Vec %i sphereHarm[%u](m=%i) = %f + %fi\n",i,P2-sphereharm_.begin(),midx,Dcomplex[0],Dcomplex[1]);
      sphereharm_.push_back(Dcomplex[0]);
      sphereharm_.push_back(Dcomplex[1]);
    }
    VXYZ += 6;
  }
}

// DataSet_Vector::FillData()
int DataSet_Vector::FillData(ComplexArray& dest, int midx) {
  if ( abs(midx) > order_ ) {
    mprinterr("Internal Error: Vector %s: FillData called with m=%i (max order= %i)\n",
              Legend().c_str(), midx, order_);
    return 1;
  }
  int p2blocksize = 2 * (2 * order_ + 1);
  double *CF = dest.CAptr();
  for ( std::vector<double>::iterator sidx = sphereharm_.begin() + 2 * (midx + order_);
                                      sidx < sphereharm_.end();
                                      sidx += p2blocksize)
  {
    *(CF++) = *sidx;
    *(CF++) = *(sidx+1);
  }
  return 0;
}
