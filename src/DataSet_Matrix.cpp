#include "DataSet_Matrix.h"

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
  type_(MATRIX_NULL)
{}

// DESTRUCTOR
DataSet_Matrix::~DataSet_Matrix() {
  if (mat_!=0) delete[] mat_;
  if (vect_!=0) delete[] vect_;
  if (vect2_!=0) delete[] vect2_;
  if (mass_!=0) delete[] mass_;
}

