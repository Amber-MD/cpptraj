// DataSet_Matrix
#include "DataSet_Matrix.h"
#include <cstring> // memcpy

using namespace std;

// CONSTRUCTOR
DataSet_Matrix::DataSet_Matrix() { 
  n_rows = 0;
  m_cols = 0;
  msize = 0;
  MSIZE_BYTES = 0;
  dType = MATRIX;
}

// CONSTRUCTOR
/// Take number of rows and columns.
DataSet_Matrix::DataSet_Matrix(int n, int m) {
  n_rows = n;
  m_cols = m;
  msize = n_rows * m_cols;
  MSIZE_BYTES = msize * sizeof(double);
  dType = MATRIX;
}

// DESTRUCTOR
DataSet_Matrix::~DataSet_Matrix() {
  for (map<int,double*>::iterator mat = matrices.begin();
                                       mat!= matrices.end();
                                       mat++)
    delete[] (*mat).second;
}

// COPY CONSTRUCTOR
DataSet_Matrix::DataSet_Matrix(const DataSet_Matrix &rhs) {
  double *mIn;
  map<int,double*>::iterator mposition = matrices.begin();
 
  n_rows = rhs.n_rows;
  m_cols = rhs.m_cols;
  msize = rhs.msize;
  dType = MATRIX;
  MSIZE_BYTES = msize * sizeof(double);
  if (msize > 0) {
    for (map<int,double*>::const_iterator mat = rhs.matrices.begin();
                                               mat!= rhs.matrices.end();
                                               mat++)
    {
      mIn = new double[ rhs.msize ];
      memcpy(mIn, (*mat).second, MSIZE_BYTES);
      mposition = matrices.insert( mposition, pair<int,double*>((*mat).first, mIn) );
    }
  }
}

// DataSet_Matrix::operator=()
DataSet_Matrix &DataSet_Matrix::operator=(const DataSet_Matrix &rhs) {
  double *mIn;
  map<int,double*>::iterator mposition;
  // Check for self assignment
  if (this == &rhs) return *this;
  // Deallocate
  for (map<int,double*>::iterator mat = matrices.begin();
                                       mat!= matrices.end();
                                       mat++)
    delete[] (*mat).second;
  matrices.clear();
  mposition = matrices.begin();
  // Allocate and copy
  n_rows = rhs.n_rows;
  m_cols = rhs.m_cols;
  msize = rhs.msize;
  MSIZE_BYTES = msize * sizeof(double);
  if (msize > 0) {
    for (map<int,double*>::const_iterator mat = rhs.matrices.begin();
                                               mat!= rhs.matrices.end();
                                               mat++)
    {
      mIn = new double[ rhs.msize ];
      memcpy(mIn, (*mat).second, MSIZE_BYTES);
      mposition = matrices.insert( mposition, pair<int,double*>((*mat).first, mIn) );
    }
  }
  // return
  return *this;
}

// DataSet_Matrix::Add()
void DataSet_Matrix::Add(int idx, void *void_matrixIn) {
  double *matrixIn = (double*) void_matrixIn;
  // It is assumed that matrix is of size msize
  double *matrix = new double[ msize ];
  memcpy(matrix, matrixIn, MSIZE_BYTES);

  matrices[idx] = matrix;
}

// DataSet_Matrix::Get()
int DataSet_Matrix::Get(void *void_matrixOut, int idx) {
  if (void_matrixOut==NULL) return 1;
  double *matrixOut = (double*) void_matrixOut;
  // It is assumed that matrix is of size msize
  map<int,double*>::iterator mposition = matrices.find( idx );
  if (mposition == matrices.end()) return 1;
  memcpy(matrixOut, (*mposition).second, MSIZE_BYTES);
  return 0;
}
 
