// TriangleMatrix
#include "TriangleMatrix.h"
#include <cfloat>
#include <cstdlib>
#include "CpptrajStdio.h"

// CONSTRUCTOR
TriangleMatrix::TriangleMatrix() {
  elements=NULL;
  nrows=0;
  nelements=0;
  currentElement=0;
  ignore=NULL;
}

// DESTRUCTOR
TriangleMatrix::~TriangleMatrix() {
  if (elements!=NULL) free(elements);
  if (ignore!=NULL) free(ignore);
}

/* TriangleMatrix::Setup()
 * Set matrix up based on the given size of 1 side of the square matrix.
 * Set the current element to 0.
 */
int TriangleMatrix::Setup(int sizeIn) {
  nrows = sizeIn;
  // Use half square matrix minus the diagonal
  nelements = ( (nrows * nrows) - nrows ) / 2;
  elements = (double*) realloc(elements, nelements * sizeof(double));
  ignore = (bool*) realloc(ignore, nrows * sizeof(bool));
  for (int n=0; n<nrows; n++) ignore[n]=false;
  currentElement=0;
  if (elements==NULL) return 1;
  return 0;
}

/* TriangleMatrix::Ignore()
 * Indicate given row/col should be ignored.
 */
void TriangleMatrix::Ignore(int row) {
  ignore[row] = true;
}

/* TriangleMatrix::Copy()
 */
TriangleMatrix *TriangleMatrix::Copy() {
  TriangleMatrix *copy;
  copy = new TriangleMatrix();
  copy->Setup( this->nrows );
  for (int N = 0; N < this->nelements; N++)
    copy->AddElement( this->elements[N] );
  for (int N = 0; N < this->nrows; N++)
    if (this->ignore[N]) copy->Ignore(N);
  return copy;
}

/* TriangleMatrix::AddElement()
 * Add the input value to the element array and increment currentElement.
 * Return 1 on success, 0 if no more elements can be added.
 */
int TriangleMatrix::AddElement(double elementIn) {
  if (currentElement>=nelements) return 0;
  elements[currentElement] = elementIn;
  currentElement++;
  return 1;
}

/* TriangleMatrix::calcIndex()
 * Calculate index in elements array for given row and column.
 */
int TriangleMatrix::calcIndex(int iIn, int jIn) {
  int i, j, i1;

  if (iIn > jIn) {
    j = iIn;
    i = jIn;
  } else {
    i = iIn;
    j = jIn;
  } 
 
  i1 = i + 1;
  return ( ( (nrows * i) - ((i1 * i) / 2) ) + j - i1 );
}

/* TriangleMatrix::AddElement()
 * Set element at specified row and column.
 */
void TriangleMatrix::SetElement(int iIn, int jIn, double elementIn) {
  int idx;

  idx = calcIndex(iIn, jIn);

  elements[idx] = elementIn;
}

/* TriangleMatrix::GetElement()
 * Get the element at specified row and column.
 */
double TriangleMatrix::GetElement(int iIn, int jIn) {
  int idx;
  
  idx = calcIndex(iIn, jIn);

  return elements[idx];
}

/* TriangleMatrix::FindMin()
 * Find the minimum, set row and column. 
 */
double TriangleMatrix::FindMin(int *iOut, int *jOut) {
  double min;
  int iVal, jVal;

  *iOut = -1;
  *jOut = -1;
  if (elements==NULL || nelements < 1) return 0.0;

  iVal = 0;
  jVal = 1;
  min = DBL_MAX;
  for (int idx = 0; idx < nelements; idx++) {
    // If we dont care about this row/col, just increment
    if (ignore[iVal] || ignore[jVal]) {
      // DEBUG
      //mprintf("\t\tIgnoring %i %i\n",iVal,jVal);
      // Increment indices
      jVal++;
      if (jVal == nrows) {
        iVal++;
        jVal = iVal + 1;
      }
    // Otherwise search for minimum
    } else {
      if ( elements[idx] < min ) {
        min = elements[idx];
        *iOut = iVal;
        *jOut = jVal;
      }
      // Increment indices
      jVal++;
      if (jVal == nrows) {
        iVal++;
        jVal = iVal + 1;
      }
    }
  }
  return min;
}

/* TriangleMatrix::PrintElements()
 */
void TriangleMatrix::PrintElements() {
  int iVal = 0;
  int jVal = 1;

  for (int idx = 0; idx < nelements; idx++) {
    if (!ignore[iVal] && !ignore[jVal])
      mprintf("\t%i %i %8.3lf\n",iVal,jVal,elements[idx]);
    // Increment indices
    jVal++;
    if (jVal == nrows) {
      iVal++;
      jVal = iVal + 1;
    }
  }
}
