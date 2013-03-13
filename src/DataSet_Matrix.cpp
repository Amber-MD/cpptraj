#include <cstring> // memset
#include "DataSet_Matrix.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet_Matrix::DataSet_Matrix() :
  DataSet(MATRIX, 6, 3, 2),
  mat_(0),
  vect_(0),
  vect2_(0),
  mass_(0),
  type_(NO_OP),
  matsize_(0),
  nrows_(0),
  ncols_(0),
  nelts_(0),
  vectsize_(0),
  calcIndex(0)
{}

// DESTRUCTOR
DataSet_Matrix::~DataSet_Matrix() {
  if (mat_!=0) delete[] mat_;
  if (vect_!=0) delete[] vect_;
  if (vect2_!=0) delete[] vect2_;
  if (mass_!=0) delete[] mass_;
}

/// Strings corresponding to MatrixType enumerated type.
const char* DataSet_Matrix::MatrixTypeString[] = {
  "UNDEFINED",
  "distance matrix",
  "covar matrix",
  "mass weighted covar matrix",
  "correlation matrix",
  "distance covar matrix",
  "idea matrix",
  "ired matrix"
};

const char* DataSet_Matrix::MatrixOutputString[] = {
  "UNKNOWN",
  "DIST",
  "COVAR",
  "MWCOVAR",
  "CORREL",
  "DISTCOVAR",
  "IDEA",
  "IRED"
};

/** Determine matrix type from arglist. Default to DIST. */
DataSet_Matrix::MatrixType DataSet_Matrix::TypeFromArg(ArgList& actionArgs) {
  MatrixType mtype = DIST;
  if (actionArgs.hasKey("distcovar"))
    mtype = DISTCOVAR;
  else if (actionArgs.hasKey("mwcovar"))
    mtype = MWCOVAR;
  else if (actionArgs.hasKey("dist"))
    mtype = DIST;
  else if (actionArgs.hasKey("covar"))
    mtype = COVAR;
  else if (actionArgs.hasKey("correl"))
    mtype = CORREL;
  else if (actionArgs.hasKey("idea"))
    mtype = IDEA;
  else if (actionArgs.hasKey("ired"))
    mtype = IRED;
  return mtype;
}

// DataSet_Matrix::AllocateVectors()
int DataSet_Matrix::AllocateVectors(size_t sizeIn) {
  if (vect_==0) {
    vectsize_ = sizeIn;
    if (vectsize_ > 0) {
      vect_ = new double[ vectsize_ ];
      memset(vect_, 0, vectsize_*sizeof(double));
      vect2_ = new double[ vectsize_ ];
      memset(vect2_, 0, vectsize_*sizeof(double));
    }
  } else {
    // Already allocated. Make sure size does not change.
    if (sizeIn != vectsize_) {
      mprinterr("Error: DataSet_Matrix: Attempting to reallocate matrix vector with \n");
      mprinterr("Error: different size. Original size= %u, new size= %u\n", vectsize_, sizeIn);
      mprinterr("Error: This can occur when different #s of atoms are selected in\n");
      mprinterr("Error: different topology files.\n");
      return 1;
    }
  }
  return 0;
}

// DataSet_Matrix::AllocateMatrix()
int DataSet_Matrix::AllocateMatrix(size_t ncolsIn, size_t nrowsIn, size_t neltsIn, size_t sizeIn)
{
  if (mat_ == 0) {
    matsize_ = sizeIn;
    if (matsize_ > 0) {
      mat_ = new double[ matsize_ ];
      memset(mat_, 0, matsize_ * sizeof(double));
    }
    ncols_ = ncolsIn;
    nrows_ = nrowsIn;
    nelts_ = neltsIn;
    if (nrows_ == 0) 
      calcIndex = calcHalfIndex;
    else
      calcIndex = calcFullIndex;
  } else {
    // Already allocated. Make sure size does not change.
    if (matsize_ != sizeIn) {
      mprinterr("Error: DataSet_Matrix: Attempting to reallocate matrix with different size.\n");
      mprinterr("Error: Original size= %u, new size= %u\n", matsize_, sizeIn);
      mprinterr("Error: This can occur when different #s of atoms are selected in\n");
      mprinterr("Error: different topology files.\n");
      return 1;
    }
  }
  return 0;
}

// DataSet_Matrix::StoreMass()
void DataSet_Matrix::StoreMass(std::vector<double> const& massIn) {
  // Store masses for mass-weighted covariance
  if (mass_ == 0) {
    mass_ = new double[ massIn.size() ];
    double* Mass = mass_;
    for (std::vector<double>::const_iterator mass = massIn.begin();
                                             mass != massIn.end(); ++mass)
      *(Mass++) = *mass;
  }
}

// DataSet_Matrix::calcFullIndex()
size_t DataSet_Matrix::calcFullIndex(size_t ncols, size_t i, size_t j) { 
  return ( (i*ncols)+j ); 
}

// DataSet_Matrix::calcHalfIndex()
size_t DataSet_Matrix::calcHalfIndex(size_t ncols, size_t row, size_t col) {
  size_t i, j;
  if (row>col) {
    i = col;
    j = row;
  } else {
    i = row;
    j = col;
  }
  //mprintf("CDBG:\ti=%i j=%i N=%i idx=%i\n",i,j,Nelt_,
  //        (i * Nelt_ - (i * (i-1) / 2) + (j - i)));
  return (i * ncols - (i * (i-1) / 2) + (j - i));
}

// DataSet_Matrix::Write2D()
void DataSet_Matrix::Write2D( CpptrajFile& outfile, int x, int y ) {
  size_t index = calcIndex(ncols_, x, y);
  if (index >= matsize_)
    outfile.Printf(data_format_, 0.0);
  else
    outfile.Printf(data_format_, mat_[index]);
}

// DataSet_Matrix::GetDimensions()
void DataSet_Matrix::GetDimensions( std::vector<int>& vIn ) {
  vIn.resize( 2 );
  if (nrows_ > 0) {
    vIn[0] = (int)ncols_;
    vIn[1] = (int)nrows_;
  } else {
    vIn[0] = (int)ncols_;
    vIn[1] = (int)ncols_;
  }
}

// DataSet_Matrix::DivideBy()
void DataSet_Matrix::DivideBy(double dval) {
  if (vect_!=0) {
    for (size_t i = 0; i < vectsize_; ++i) {
      vect_[i] /= dval;
      vect2_[i] /= dval;
    }
  }
  for (size_t i = 0; i < matsize_; ++i)
    mat_[i] /= dval;
}

/** Calc <riri> - <ri><ri> */
void DataSet_Matrix::Vect2MinusVect() {
  for (size_t i = 0; i < vectsize_; ++i) {
    //mprintf("CDBG:\tvect2[%i] = %lf\n",i,vect2_[i]);
    vect2_[i] -= (vect_[i]*vect_[i]);
    //mprintf("CDBG:\tvect2[%i] = %lf\n",i,vect2_[i]);
  }
}

