// TriangleMatrix
#include "TriangleMatrix.h"
#include <cfloat>
#include "CpptrajStdio.h"
#include <cstring> //memcpy
#include <cstdio> // save and load

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
  if (elements!=NULL) delete[] elements;
  if (ignore!=NULL) delete[] ignore;
}

/* TriangleMatrix::SaveFile()
 */
int TriangleMatrix::SaveFile(char *filename) {
  FILE *outfile;
  char magic[4];

  magic[0]='C'; // Cpptraj
  magic[1]='T'; // Triangle
  magic[2]='M'; // Matrix
  magic[3]=0;   // Version

  outfile = fopen(filename,"wb");
  if (outfile==NULL) {
    mprinterr("Error: TriangleMatrix::SaveFile: Could not open file %s\n",filename);
    return 1;
  }
  // Write magic byte
  fwrite(magic, sizeof(char), 4, outfile);
  // Write nrows
  fwrite(&nrows, sizeof(int), 1, outfile);
  // Write number of elements
  fwrite(&nelements, sizeof(int), 1, outfile);
  // Write elements
  fwrite(elements, sizeof(float), nelements, outfile);
  // Write ignore
  //fwrite(ignore, sizeof(bool), nrows, outfile);

  fclose(outfile);
  return 0;
}

/* TriangleMatrix::LoadFile
 */
int TriangleMatrix::LoadFile(char *filename, int sizeIn) {
  FILE *infile;
  char magic[4];

  infile = fopen(filename, "rb");
  if (infile==NULL) {
    mprinterr("Error: TriangleMatrix::LoadFile: Could not open file %s\n",filename);
    return 1;
  }
  // Read and check magic byte
  fread(magic, sizeof(char), 4, infile);
  if (magic[0] != 'C' ||
      magic[1] != 'T' ||
      magic[2] != 'M' ||
      magic[3] != 0     ) {
    mprinterr("Error: TriangleMatrix::LoadFile: File %s magic number not CTM0!\n",filename);
    fclose(infile);
    return 1;
  }
  // Read nrows
  fread(&nrows, sizeof(int), 1, infile);
  // If number of rows is not what was expected, abort
  if (nrows != sizeIn) {
    mprintf("Warning: TriangleMatrix::LoadFile: File %s has %i rows, expected %i.\n",nrows,sizeIn);
    fclose(infile);
    return 1;
  }
  // Read number of elements
  fread(&nelements, sizeof(int), 1, infile);
  if (elements!=NULL) delete[] elements;
  elements = new float[ nelements ];
  // Read elements
  fread(elements,sizeof(float), nelements, infile);
  // Setup ignore array
  if (ignore!=NULL) delete[] ignore;
  ignore = new bool[ nrows ];
  for (int n=0; n<nrows; n++) ignore[n]=false;
  currentElement=0;

  fclose(infile);
  return 0;
}

/* TriangleMatrix::Setup()
 * Set matrix up based on the given size of 1 side of the square matrix.
 * Set the current element to 0.
 */
int TriangleMatrix::Setup(int sizeIn) {
  nrows = sizeIn;
  // Use half square matrix minus the diagonal
  nelements = ( (nrows * nrows) - nrows ) / 2;
  if (elements!=NULL) delete[] elements;
  elements = new float[ nelements ];
  // Setup ignore array
  if (ignore!=NULL) delete[] ignore;
  ignore = new bool[ nrows ];
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

/* TriangleMatrix::operator=()
 */
TriangleMatrix &TriangleMatrix::operator=(const TriangleMatrix &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;

  // Deallocate
  if (elements!=NULL) delete[] elements;
  if (ignore!=NULL) delete[] ignore;

  // Allocate
  nelements = rhs.nelements;
  nrows = rhs.nrows;
  currentElement = rhs.currentElement;
  elements = new float[ nelements ];
  ignore = new bool[ nrows ];

  // Copy
  memcpy(elements, rhs.elements, nelements * sizeof(float));
  memcpy(ignore,   rhs.ignore,   nrows * sizeof(bool));

  // Return *this
  return *this;
}

/* TriangleMatrix::AddElement()
 * Add the input value to the element array and increment currentElement.
 * Return 1 on success, 0 if no more elements can be added.
 */
int TriangleMatrix::AddElement(double elementIn) {
  if (currentElement>=nelements) return 0;
  elements[currentElement] = (float) elementIn;
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

  elements[idx] = (float) elementIn;
}

/* TriangleMatrix::GetElement()
 * Get the element at specified row and column.
 */
double TriangleMatrix::GetElement(int iIn, int jIn) {
  int idx;
  
  idx = calcIndex(iIn, jIn);

  return (double)elements[idx];
}

/* TriangleMatrix::FindMin()
 * Find the minimum, set row and column. 
 */
double TriangleMatrix::FindMin(int *iOut, int *jOut) {
  float min;
  int iVal, jVal;

  *iOut = -1;
  *jOut = -1;
  if (elements==NULL || nelements < 1) return 0.0;

  iVal = 0;
  jVal = 1;
  min = FLT_MAX;
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
  return (double)min;
}

/* TriangleMatrix::PrintElements()
 */
void TriangleMatrix::PrintElements() {
  int iVal = 0;
  int jVal = 1;

  for (int idx = 0; idx < nelements; idx++) {
    if (!ignore[iVal] && !ignore[jVal])
      mprintf("\t%i %i %8.3f\n",iVal,jVal,elements[idx]);
    // Increment indices
    jVal++;
    if (jVal == nrows) {
      iVal++;
      jVal = iVal + 1;
    }
  }
}
