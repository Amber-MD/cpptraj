// TriangleMatrix
#include <cstring> //memcpy
#include "TriangleMatrix.h"

// CONSTRUCTOR
TriangleMatrix::TriangleMatrix() :
  DataSet(TRIMATRIX, 12, 4, 2),
  elements_(0),
  nrows_(0),
  nelements_(0),
  currentElement_(0)
{}

// CONSTRUCTOR
TriangleMatrix::TriangleMatrix(int sizeIn) :
  DataSet(TRIMATRIX, 12, 4, 2),
  elements_(0),
  nrows_(sizeIn),
  nelements_( (size_t)( (nrows_*(nrows_ - 1L)) / 2L) ),
  currentElement_(0)
{
  if (nrows_ > 0)
    elements_ = new float[ nelements_ ];
}

// DESTRUCTOR
TriangleMatrix::~TriangleMatrix() {
  if (elements_!=0) delete[] elements_;
}

// COPY CONSTRUCTOR
TriangleMatrix::TriangleMatrix(const TriangleMatrix &rhs) :
  DataSet( rhs )
{
  nelements_ = rhs.nelements_;
  nrows_ = rhs.nrows_;
  currentElement_ = rhs.currentElement_;
  if (nrows_ > 0) {
    elements_ = new float[ nelements_ ];
    memcpy(elements_, rhs.elements_, nelements_ * sizeof(float));
  } else
    elements_ = 0;
}

// TriangleMatrix::operator=()
TriangleMatrix &TriangleMatrix::operator=(const TriangleMatrix &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  DataSet::operator=(rhs);
  // Deallocate
  if (elements_!=0) delete[] elements_;
  // Allocate
  nelements_ = rhs.nelements_;
  nrows_ = rhs.nrows_;
  currentElement_ = rhs.currentElement_;
  if (nrows_ > 0) {
    elements_ = new float[ nelements_ ];
    // Copy
    memcpy(elements_, rhs.elements_, nelements_ * sizeof(float));
  } else
    elements_ = 0;
  return *this;
}

// TriangleMatrix::Setup()
/** Set matrix up based on the given size of 1 side of the square matrix.
  * Set the current element to 0.
  */
int TriangleMatrix::Setup(int sizeIn) {
  nrows_ = sizeIn;
  size_t ROWS = (size_t) nrows_;
  // Use half square matrix minus the diagonal
  nelements_ = ( ROWS * (ROWS - 1) ) / 2; 
  if (elements_!=0) delete[] elements_;
  if (nelements_ > 0)
    elements_ = new float[ nelements_ ];
  else
    elements_ = 0;
  currentElement_ = 0;
  if (elements_==0) return 1;
  return 0;
}

// TriangleMatrix::AddElement()
/** Add the input double to the element array and increment currentElement.
  * \return 1 on success, 0 if no more elements can be added.
  */
int TriangleMatrix::AddElement(double elementIn) {
  if (currentElement_>=nelements_) return 0;
  elements_[currentElement_] = (float) elementIn;
  ++currentElement_;
  return 1;
}

// TriangleMatrix::AddElement()
/** Add the input float to the element array and increment currentElement.
  * \return 1 on success, 0 if no more elements can be added.
  */
int TriangleMatrix::AddElement(float elementIn) {
  if (currentElement_>=nelements_) return 0;
  elements_[currentElement_] = elementIn;
  ++currentElement_;
  return 1;
}

// TriangleMatrix::calcIndex()
/** Calculate index in elements array for given row and column.
  * SHOULD NEVER BE CALLED WITH iIn == jIn!
  */
int TriangleMatrix::calcIndex(int iIn, int jIn) const {
  int i, j;
  if (iIn > jIn) {
    j = iIn;
    i = jIn;
  } else {
    i = iIn;
    j = jIn;
  } 
  int i1 = i + 1;
  return ( ( (nrows_ * i) - ((i1 * i) / 2) ) + j - i1 );
}

// TriangleMatrix::SetElement()
/** Set element at specified row and column. */
void TriangleMatrix::SetElement(int iIn, int jIn, double elementIn) {
  if (iIn == jIn) return;
  int idx = calcIndex(iIn, jIn);
  elements_[idx] = (float) elementIn;
}

// TriangleMatrix::SetElementF()
/** Set element at specified row and column. */
void TriangleMatrix::SetElementF(int iIn, int jIn, float elementIn) {
  if (iIn == jIn) return;
  int idx = calcIndex(iIn, jIn);
  elements_[idx] = elementIn;
}


// TriangleMatrix::GetElement()
/** Get the element at specified row and column as a double.
  */
double TriangleMatrix::GetElement(int iIn, int jIn) const {
  if (iIn == jIn) return 0;
  int idx = calcIndex(iIn, jIn);
  return (double)elements_[idx];
}

// TriangleMatrix::GetElementF()
/** Get the element at specified row and column. */
float TriangleMatrix::GetElementF(int iIn, int jIn) const {
  if (iIn == jIn) return 0;
  int idx = calcIndex(iIn, jIn);
  return elements_[idx];
}

void TriangleMatrix::Write2D( CpptrajFile& outfile, int x, int y ) {
  if ( x==y || x < 0 || y < 0 || x >= nrows_ || y >= nrows_ ) 
    outfile.Printf(data_format_, 0.0);
  else {
    int index = calcIndex(x, y);
    outfile.Printf(data_format_, elements_[index]);
  }
}

void TriangleMatrix::GetDimensions( std::vector<int>& vIn ) {
  vIn.assign( 2, nrows_ );
}
