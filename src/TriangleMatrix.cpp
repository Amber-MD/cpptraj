// TriangleMatrix
#include <cfloat> // FLT_MAX
#include <cstring> //memcpy
#include <cstdio> // save and load
#include "TriangleMatrix.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
TriangleMatrix::TriangleMatrix() : 
  elements_(0),
  nrows_(0),
  nelements_(0),
  currentElement_(0),
  ignore_(0)
{
  // DataSet-specific vars
  width_ = 12;
  precision_ = 4;
  dType_ = TRIMATRIX;
  SetDataSetFormat(false);
  dim_ = 2;
}

// CONSTRUCTOR
TriangleMatrix::TriangleMatrix(int sizeIn) :
  elements_(0),
  nrows_(sizeIn),
  nelements_( (size_t)( (nrows_*(nrows_ - 1)) / 2) ),
  currentElement_(0),
  ignore_(0)
{
  elements_ = new float[ nelements_ ];
  ignore_ = new bool[ nrows_ ];
  for (int n=0; n<nrows_; n++) 
    ignore_[n]=false;
  // DataSet-specific vars
  width_ = 12;
  precision_ = 4;
  dType_ = TRIMATRIX;
  SetDataSetFormat(false);
  dim_ = 2;
}

// DESTRUCTOR
TriangleMatrix::~TriangleMatrix() {
  if (elements_!=0) delete[] elements_;
  if (ignore_!=0) delete[] ignore_;
}

// COPY CONSTRUCTOR
TriangleMatrix::TriangleMatrix(const TriangleMatrix &rhs) :
  DataSet( rhs )
{
  nelements_ = rhs.nelements_;
  nrows_ = rhs.nrows_;
  currentElement_ = rhs.currentElement_;
  elements_ = new float[ nelements_ ];
  ignore_ = new bool[ nrows_ ];
  memcpy(elements_, rhs.elements_, nelements_ * sizeof(float));
  memcpy(ignore_,   rhs.ignore_,   nrows_ * sizeof(bool));
}

// TriangleMatrix::operator=()
TriangleMatrix &TriangleMatrix::operator=(const TriangleMatrix &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  DataSet::operator=(rhs);
  // Deallocate
  if (elements_!=0) delete[] elements_;
  if (ignore_!=0) delete[] ignore_;

  // Allocate
  nelements_ = rhs.nelements_;
  nrows_ = rhs.nrows_;
  currentElement_ = rhs.currentElement_;
  elements_ = new float[ nelements_ ];
  ignore_ = new bool[ nrows_ ];

  // Copy
  memcpy(elements_, rhs.elements_, nelements_ * sizeof(float));
  memcpy(ignore_,   rhs.ignore_,   nrows_ * sizeof(bool));

  // Return *this
  return *this;
}

// TriangleMatrix::SaveFile()
/** Save the matrix to a binary file. Format is 
  *   Data: [4*char][int][int][nelements*float]
  *   Vars: ['C''T''M'0][nrows][nelements][elements]
  */
int TriangleMatrix::SaveFile(const char *filename) {
  char magic[4];

  magic[0]='C'; // Cpptraj
  magic[1]='T'; // Triangle
  magic[2]='M'; // Matrix
  magic[3]=0;   // Version

  FILE* outfile = fopen(filename,"wb");
  if (outfile==NULL) {
    mprinterr("Error: TriangleMatrix::SaveFile: Could not open file %s\n",filename);
    return 1;
  }
  // Write magic byte
  fwrite(magic, sizeof(char), 4, outfile);
  // Write nrows
  fwrite(&nrows_, sizeof(int), 1, outfile);
  // Write number of elements
  // NOTE: Make long int?
  int Ntemp = (int) nelements_;
  fwrite(&Ntemp, sizeof(int), 1, outfile);
  // Write elements
  fwrite(elements_, sizeof(float), nelements_, outfile);
  // Write ignore
  //fwrite(ignore_, sizeof(bool), nrows_, outfile);

  fclose(outfile);
  return 0;
}

// TriangleMatrix::LoadFile
int TriangleMatrix::LoadFile(const char *filename, int sizeIn) {
  char magic[4];

  FILE* infile = fopen(filename, "rb");
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
    mprinterr("Error: TriangleMatrix::LoadFile: File %s magic number not CTM0 [%c%c%c%c]!\n",
              filename, magic[0], magic[1], magic[2], magic[3]);
    fclose(infile);
    return 1;
  }
  // Read nrows
  fread(&nrows_, sizeof(int), 1, infile);
  // If number of rows is not what was expected, abort
  if (nrows_ != sizeIn) {
    mprinterr("Error: TriangleMatrix::LoadFile: File %s has %i rows, expected %i.\n",
            filename, nrows_, sizeIn);
    fclose(infile);
    return 1;
  }
  // Read number of elements
  // NOTE: Read long int?
  int Ntemp = 0;
  fread(&Ntemp, sizeof(int), 1, infile);
  nelements_ = (size_t) Ntemp;
  if (elements_!=0) delete[] elements_;
  elements_ = new float[ nelements_ ];
  // Read elements
  fread(elements_,sizeof(float), nelements_, infile);
  // Setup ignore array
  if (ignore_!=0) delete[] ignore_;
  ignore_ = new bool[ nrows_ ];
  for (int n=0; n<nrows_; n++) ignore_[n]=false;
  currentElement_=0;

  fclose(infile);
  return 0;
}

// TriangleMatrix::Setup()
/** Set matrix up based on the given size of 1 side of the square matrix.
  * Set the current element to 0.
  */
int TriangleMatrix::Setup(int sizeIn) {
  nrows_ = sizeIn;
  size_t ROWS = (size_t) nrows_;
  // Use half square matrix minus the diagonal
  //nelements = ( (nrows * nrows) - nrows ) / 2;
  nelements_ = ( ROWS * (ROWS - 1) ) / 2; 
  if (elements_!=0) delete[] elements_;
  //mprintf("DEBUG: TriangleMatrix::Setup(%i) nrows=%i nelements=%lu\n",sizeIn,nrows,nelements);
  elements_ = new float[ nelements_ ];
  // Setup ignore array
  if (ignore_!=0) delete[] ignore_;
  ignore_ = new bool[ nrows_ ];
  for (int n=0; n<nrows_; n++) ignore_[n]=false;
  currentElement_=0;
  if (elements_==0) return 1;
  return 0;
}

// TriangleMatrix::Ignore()
/** Indicate given row/col should be ignored. */
void TriangleMatrix::Ignore(int row) {
  ignore_[row] = true;
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
  return ( ( (nrows_ * i) - ((i1 * i) / 2) ) + j - i1 );
}

// TriangleMatrix::SetElement()
/** Set element at specified row and column. */
void TriangleMatrix::SetElement(int iIn, int jIn, double elementIn) {
  int idx;

  if (iIn == jIn) return;

  idx = calcIndex(iIn, jIn);

  elements_[idx] = (float) elementIn;
}

// TriangleMatrix::SetElementF()
/** Set element at specified row and column. */
void TriangleMatrix::SetElementF(int iIn, int jIn, float elementIn) {
  int idx;

  if (iIn == jIn) return;

  idx = calcIndex(iIn, jIn);

  elements_[idx] = elementIn;
}


// TriangleMatrix::GetElement()
/** Get the element at specified row and column as a double.
  */
double TriangleMatrix::GetElement(int iIn, int jIn) {
  int idx;
 
  if (iIn == jIn) return 0;
 
  idx = calcIndex(iIn, jIn);

  return (double)elements_[idx];
}

// TriangleMatrix::GetElementF()
/** Get the element at specified row and column. */
float TriangleMatrix::GetElementF(int iIn, int jIn) {
  int idx;
  
  if (iIn == jIn) return 0;
 
  idx = calcIndex(iIn, jIn);

  return elements_[idx];
}

// TriangleMatrix::FindMin()
/** Find the minimum, set row and column. 
  */
double TriangleMatrix::FindMin(int *iOut, int *jOut) {
  float min;
  int iVal, jVal;

  *iOut = -1;
  *jOut = -1;
  if (elements_==0 || nelements_ < 1) return 0.0;

  iVal = 0;
  jVal = 1;
  min = FLT_MAX;
  for (size_t idx = 0; idx < nelements_; idx++) {
    // If we dont care about this row/col, just increment
    if (ignore_[iVal] || ignore_[jVal]) {
      // DEBUG
      //mprintf("\t\tIgnoring %i %i\n",iVal,jVal);
      // Increment indices
      jVal++;
      if (jVal == nrows_) {
        iVal++;
        jVal = iVal + 1;
      }
    // Otherwise search for minimum
    } else {
      if ( elements_[idx] < min ) {
        min = elements_[idx];
        *iOut = iVal;
        *jOut = jVal;
      }
      // Increment indices
      jVal++;
      if (jVal == nrows_) {
        iVal++;
        jVal = iVal + 1;
      }
    }
  }
  return (double)min;
}

// TriangleMatrix::PrintElements()
void TriangleMatrix::PrintElements() {
  int iVal = 0;
  int jVal = 1;

  for (size_t idx = 0; idx < nelements_; idx++) {
    if (!ignore_[iVal] && !ignore_[jVal])
      mprintf("\t%i %i %8.3f\n",iVal,jVal,elements_[idx]);
    // Increment indices
    jVal++;
    if (jVal == nrows_) {
      ++iVal;
      jVal = iVal + 1;
    }
  }
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
