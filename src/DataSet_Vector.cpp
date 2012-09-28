#include <cstring> // memcpy
#include <cmath> // sqrt
#include "DataSet_Vector.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_Vector::DataSet_Vector() :
  DataSet(VECTOR, 8, 4),
  totalidx_(0),
  currentidx_(0),
  xyz_(0),
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
    memset(newxyz+totalidx_, 0,      3000 * sizeof(double));
    delete[] xyz_;
  }
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
    mprinterr("Eerror: DataSet_Vector: Frame %i is out of range.\n",frameIn);
    return;
  }
  cbuffer.Printf(data_format_, xyz_[idx  ]);
  cbuffer.Printf(data_format_, xyz_[idx+1]);
  cbuffer.Printf(data_format_, xyz_[idx+2]);
  cbuffer.Printf(data_format_, xyz_[idx+3]);
  cbuffer.Printf(data_format_, xyz_[idx+4]);
  cbuffer.Printf(data_format_, xyz_[idx+5]);
  cbuffer.Printf(data_format_, xyz_[idx  ] + xyz_[idx+3]);
  cbuffer.Printf(data_format_, xyz_[idx+1] + xyz_[idx+4]);
  cbuffer.Printf(data_format_, xyz_[idx+2] + xyz_[idx+5]);
}

// VectorType::sphericalHarmonics()
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
    } else if (std::abs(m) == 1) {
      D[0] = -m * SH21 * x * z * ri * ri;
      D[1] =     -SH21 * y * z * ri * ri;
    } else {
      D[0] = SH22 * (x*x - y*y) * ri * ri;
      D[1] = m * SH22 * x * y * ri * ri;
    }
  }
  //mprintf("CDBG: dreal=%.10lE dimg=%.10lE\n",D[0],D[1]);
}

// DataSet_Vector::SphericalHarmonics()
double* DataSet_Vector::SphericalHarmonics(int order) {
  double Dcomplex[2];

  // Allocate space to hold complex numbers. Each frame has 2 doubles
  // (real + imaginary)  for each m = -olegendre -> +olegendre
  int n_of_vecs = currentidx_ / 6;
  int p2blocksize = 2 * ((order * 2) + 1);
  int p2size = (p2blocksize * n_of_vecs);
  double* p2cftmp = new double[ p2size ];
  memset( p2cftmp, 0, p2size * sizeof(double) );
   // Loop over all vectors, convert to spherical harmonics
  double* VXYZ = xyz_;
  double* P2 = p2cftmp;
  for (int i = 0; i < n_of_vecs; ++i) {
    // Calc vector length
    double len = sqrt(VXYZ[0]*VXYZ[0] + VXYZ[1]*VXYZ[1] + VXYZ[2]*VXYZ[2]);
    // Loop over m = -olegendre, ..., +olegendre
    for (int midx = -order; midx <= order; ++midx) {
      sphericalHarmonics(order, midx, VXYZ, len, Dcomplex);
      //mprintf("DBG: sphereHarm[%i](m=%i) = %lf + %lfi\n",i,midx,Dcomplex[0],Dcomplex[1]);
      *(P2++) = Dcomplex[0];
      *(P2++) = Dcomplex[1];
    }
    VXYZ += 3;
  }
  return p2cftmp;
}

