// TriangleMatrix
#include "TriangleMatrix.h"
#include <cstdlib>

// CONSTRUCTOR
TriangleMatrix::TriangleMatrix() {
  elements=NULL;
  Nrows=0;
  Nelements=0;
  currentElement=0;
}

// DESTRUCTOR
TriangleMatrix::~TriangleMatrix() {
  if (elements!=NULL) free(elements);
}

/* TriangleMatrix::Setup()
 * Set matrix up based on the given size of 1 side of the square matrix.
 * Set the current element to 0.
 */
int TriangleMatrix::Setup(int sizeIn) {
  Nrows = sizeIn;
  // Use half square matrix minus the diagonal
  Nelements = ( (Nrows * Nrows) - Nrows ) / 2;
  elements = (double*) realloc(elements, Nelements * sizeof(double));
  currentElement=0;
  if (elements==NULL) return 1;
  return 0;
}

/* TriangleMatrix::AddElement()
 * Add the input value to the element array and increment currentElement.
 * Return 1 on success, 0 if no more elements can be added.
 */
int TriangleMatrix::AddElement(double elementIn) {
  if (currentElement>=Nelements) return 0;
  elements[currentElement] = elementIn;
  currentElement++;
  return 1;
}

/* TriangleMatrix::GetElement()
 * Get the element at specified row and column.
 */
double TriangleMatrix::GetElement(int iIn, int jIn) {
  int i, j, i1, idx;
  
  if (iIn > jIn) {
    j = iIn;
    i = jIn;
  } else {
    i = iIn;
    j = jIn;
  } 
  
  i1 = i + 1;
  idx = ( (Nrows * i) - ((i1 * i) / 2) ) + j - i1;
  return elements[idx];
}

