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
  matsize_(0),
  nrows_(0),
  ncols_(0),
  vectsize_(0),
  snap_(0),
  type_(MATRIX_NULL),
  calcIndex(0)
{}

// DESTRUCTOR
DataSet_Matrix::~DataSet_Matrix() {
  if (mat_!=0) delete[] mat_;
  if (vect_!=0) delete[] vect_;
  if (vect2_!=0) delete[] vect2_;
  if (mass_!=0) delete[] mass_;
}

// DataSet_Matrix::MatrixAlloc()
int DataSet_Matrix::MatrixAlloc(AtomMask& mask1, AtomMask& mask2, 
                                std::vector<Atom> const& Atoms)
{
  if (mask1.None() && mask2.None()) {
    mprinterr("Error: DataSet_Matrix: MemAlloc: Masks are empty.\n");
    return 1;
  }
  int mask1tot = mask1.Nselected();
  int mask2tot = mask2.Nselected();
  if (mask1tot < mask2tot) {
    mprinterr("Error: DataSet_Matrix: # of atoms in mask1 < # of atoms in mask2\n");
    return 1;
  }
  // Allocate vector memory
  if (vect_ == 0) {
    switch( type_ ) {
      case MATRIX_DIST:      vectsize_ = 0; break;
      case MATRIX_DISTCOVAR: vectsize_ = mask1tot * (mask1tot - 1) / 2; break;
      case MATRIX_COVAR:
      case MATRIX_MWCOVAR:   vectsize_ = (mask1tot + mask2tot) * 3; break;
      default:               vectsize_ = mask1tot + mask2tot;
    }
    if (vectsize_ > 0) {
      vect_ = new double[ vectsize_ ];
      memset(vect_, 0, vectsize_*sizeof(double));
      vect2_ = new double[ vectsize_ ];
      memset(vect2_, 0, vectsize_*sizeof(double));
    }
  }
  // Determine matrix size in memory
  int matrixSize = 0;
  if (mask2tot == 0) {
    // "Upper right half" matrix, including main diagonal.
    switch( type_ ) {
      case MATRIX_DISTCOVAR:
        ncols_ = vectsize_;
        matrixSize = mask1tot * (mask1tot - 1) * (ncols_ + 1) / 4;
        break;
      case MATRIX_COVAR:
      case MATRIX_MWCOVAR: 
        ncols_ = mask1tot * 3;
        matrixSize = 3 * ncols_ * (mask1tot + 1) / 2;
        break;
      default:
        ncols_ = mask1tot;
        matrixSize = ncols_ * (mask1tot + 1) / 2;
    }
    calcIndex = calcHalfIndex;
  } else {
    // Full matrix - no MATRIX_DISTCOVAR, MATRIX_IDEA, or MATRIX_IRED possible
    switch( type_ ) {
      case MATRIX_DISTCOVAR:
      case MATRIX_IDEA:
      case MATRIX_IRED:
        return 1;
        break;
      case MATRIX_COVAR:
      case MATRIX_MWCOVAR:
        matrixSize = 9 * mask1tot * mask2tot;
        ncols_ = mask1tot * 3;
        nrows_ = mask2tot * 3;
        break;
      default:   
        matrixSize = mask1tot * mask2tot;
        ncols_ = mask1tot;
        nrows_ = mask2tot;
    }
    calcIndex = calcFullIndex;
  }
  // Allocate Matrix
  if (mat_ == 0) {
    if (matrixSize < 1) {
      mprinterr("Error: DataSet_Matrix: Size is < 1\n");
      return 1;
    }
    matsize_ = matrixSize;
    mat_ = new double[ matsize_ ];
    memset(mat_, 0, matsize_ * sizeof(double));
  } else if (matsize_ != matrixSize) {
    // Already allocated. Make sure size does not change.
    mprinterr("Error: DataSet_Matrix: Attempting to reallocate matrix with different size.\n");
    mprinterr("Error: Original size= %i, new size= %i\n", matsize_, matrixSize);
    mprinterr("Error: This can occur when different #s of atoms are selected in\n");
    mprinterr("Error: different topology files.\n");
    return 1;
  }
  // Store masses for mass-weighted covariance
  if (type_ == MATRIX_MWCOVAR && mass_ == 0) {
    mass_ = new double[ mask1.Nselected() ];
    double* Mass = mass_;
    for (AtomMask::const_iterator atomi = mask1.begin();
                                  atomi != mask1.end(); ++atomi)
      *(Mass++) = Atoms[ *atomi ].Mass();
  }
  
  return 0;
}

int DataSet_Matrix::calcFullIndex(int ncols, int i, int j) { 
  return ( (j*ncols)+i ); 
}

int DataSet_Matrix::calcHalfIndex(int ncols, int row, int col) {
  int i, j;
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

void DataSet_Matrix::Write2D( CpptrajFile& outfile, int x, int y ) {
  int index = calcIndex(ncols_, x, y);
  if (index < 0 || index >= matsize_)
    outfile.Printf(data_format_, 0.0);
  else
    outfile.Printf(data_format_, mat_[index]);
}

void DataSet_Matrix::GetDimensions( std::vector<int>& vIn ) {
  vIn.resize( 2 );
  if (nrows_ > 0) {
    vIn[0] = ncols_;
    vIn[1] = nrows_;
  } else {
    vIn[0] = ncols_;
    vIn[1] = ncols_;
  }
}

void DataSet_Matrix::AverageOverSnapshots() {
  double dsnap = (double)snap_;
  if (vect_!=0) {
    for (int i = 0; i < vectsize_; ++i) {
      vect_[i] /= dsnap;
      vect2_[i] /= dsnap;
    }
  }
  for (int i = 0; i < matsize_; ++i)
    mat_[i] /= dsnap;
}

/** Calc <riri> - <ri><ri> */
void DataSet_Matrix::Vect2MinusVect() {
  for (int i = 0; i < vectsize_*3; ++i) {
    //mprintf("CDBG:\tvect2[%i] = %lf\n",i,vect2_[i]);
    vect2_[i] -= (vect_[i]*vect_[i]);
    //mprintf("CDBG:\tvect2[%i] = %lf\n",i,vect2_[i]);
  }
}

