#include <cstring> // memcpy
#include "Matrix_2D.h"

// CONSTRUCTOR
Matrix_2D::Matrix_2D() :
  DataSet(MATRIX2D, 12, 4, 2),
  elements_(0),
  ncols_(0),
  nrows_(0),
  nelements_(0),
  currentElement_(0)
{ }

// DESTRUCTOR
Matrix_2D::~Matrix_2D() {
  if (elements_!=0) delete[] elements_;
}

Matrix_2D::Matrix_2D(const Matrix_2D& rhs) :
  elements_(0),
  ncols_( rhs.ncols_ ),
  nrows_( rhs.nrows_ ),
  nelements_( rhs.nelements_ ),
  currentElement_( rhs.currentElement_ )
{
  if (nelements_ > 0) {
    elements_ = new double[ nelements_ ];
    memcpy( elements_, rhs.elements_, nelements_*sizeof(double) );
  }
}

Matrix_2D& Matrix_2D::operator=(const Matrix_2D& rhs) {
  if (this == &rhs) return *this;
  
  if (elements_!=0) {
    delete[] elements_;
    elements_ = 0;
  }

  ncols_ = rhs.ncols_;
  nrows_ = rhs.nrows_;
  nelements_ = rhs.nelements_;
  currentElement_ = rhs.currentElement_;
  if (nelements_ > 0) {
    elements_ = new double[ nelements_ ];
    memcpy( elements_, rhs.elements_, nelements_*sizeof(double) );
  }
  return *this;
}

int Matrix_2D::Setup(int cols, int rows) {
  if (elements_!=0) {
    delete[] elements_;
    elements_ = 0;
  }
  ncols_ = cols;
  nelements_ = (size_t)ncols_;
  nrows_ = rows;
  nelements_ *= (size_t)nrows_;
  currentElement_ = 0;
  elements_ = new double[ nelements_ ];
  memset(elements_, 0, nelements_*sizeof(double) );
  return 0;
}

/** Add the input double to the element array and increment currentElement.
  * \return 1 on success, 0 if no more elements can be added.
  */
int Matrix_2D::AddElement(double elementIn) {
  if (currentElement_>=nelements_) return 0;
  elements_[currentElement_] = elementIn;
  ++currentElement_;
  return 1;
}

/** Set element at specified row and column. */
void Matrix_2D::SetElement(int iIn, int jIn, double elementIn) {
  int idx = calcIndex(iIn, jIn);

  elements_[idx] = elementIn;
}

/** Get the element at specified row and column. */
double Matrix_2D::GetElement(int iIn, int jIn) {
  int idx = calcIndex(iIn, jIn);

  return elements_[idx];
}

void Matrix_2D::Write2D( CpptrajFile& outfile, int x, int y ) {
  if ( x < 0 || y < 0 || x >= ncols_ || y >= nrows_ )
    outfile.Printf(data_format_, 0.0);
  else {
    int index = calcIndex(x, y);
    outfile.Printf(data_format_, elements_[index]);
  }
}

void Matrix_2D::GetDimensions( std::vector<int>& vIn ) {
  vIn.resize( 2 );
  vIn[0] = ncols_;
  vIn[1] = nrows_;
}
