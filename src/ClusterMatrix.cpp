#include <cfloat> // FLT_MAX
#include "ClusterMatrix.h"
#include "CpptrajStdio.h"

const unsigned char ClusterMatrix::Magic_[4] = {'C', 'T', 'M', 0};

/// CONSTRUCTOR - Set up TriangleMatrix and Ignore array.
ClusterMatrix::ClusterMatrix(int sizeIn) :
  TriangleMatrix(sizeIn),
  ignore_(sizeIn, false)
{}

// COPY CONSTRUCTOR
ClusterMatrix::ClusterMatrix(const ClusterMatrix& rhs) :
  TriangleMatrix(rhs),
  ignore_(rhs.ignore_)
{}

// ASSIGNMENT
ClusterMatrix& ClusterMatrix::operator=(const ClusterMatrix& rhs) {
  if (this == &rhs) return *this;
  TriangleMatrix::operator=(rhs);
  ignore_ = rhs.ignore_;
  return *this;
}

// ClusterMatrix::SaveFile()
/** Save the matrix to a binary file. Format is:
  *   Data: [4*char][int][int][nelements*float]
  *   Vars: ['C''T''M'0][nrows][nelements][elements]
  */
int ClusterMatrix::SaveFile(std::string const& filename) const {
  CpptrajFile outfile;
  // No stdout write allowed.
  if (filename.empty()) {
    mprinterr("Internal Error: ClusterMatrix::SaveFile called with no filename.\n");
    return 1;
  }
  if (outfile.OpenWrite(filename)) {
    mprinterr("Error: ClusterMatrix::SaveFile: Could not open %s for write.\n", filename.c_str());
    return 1;
  }
  // Write magic byte
  outfile.Write( Magic_, 4 );
  // Write nrows
  outfile.Write( &nrows_, sizeof(int) );
  // Write number of elements.
  // TODO: Make long int?
  int Ntemp = (int)nelements_;
  outfile.Write( &Ntemp, sizeof(int) );
  // Write matrix elements
  outfile.Write( elements_, Ntemp*sizeof(float) );
  // TODO: Write ignore?
  return 0;
}

// ClusterMatrix::LoadFile()
int ClusterMatrix::LoadFile(std::string const& filename, int sizeIn) {
  unsigned char magic[4];
  CpptrajFile infile;
  if (infile.OpenRead(filename)) {
    mprinterr("Error: ClusterMatrix::LoadFile: Could not open %s for read.\n", filename.c_str());
    return 1;
  }
  // Read and check magic byte
  infile.Read( magic, 4 );
  if ( magic[0]!=Magic_[0] || magic[1]!=Magic_[1] || magic[2]!=Magic_[2] ) {
    mprinterr("Error: ClusterMatrix::LoadFile: File %s is not a TriangleMatrix.\n",
              filename.c_str());
    return 1;
  }
  if (magic[3] != Magic_[3])
    mprintf("Warning: ClusterMatrix::LoadFile: %s version byte %c does not match (%c)\n",
            filename.c_str(), magic[3], Magic_[3]);
  // Read nrows
  int Ntemp = 0;
  infile.Read( &Ntemp, sizeof(int) );
  // If number of rows is not what was expected, abort
  if (Ntemp != sizeIn) {
    mprinterr("Error: ClusterMatrix::LoadFile: File %s has %i rows, expected %i.\n",
              filename.c_str(), Ntemp, sizeIn);
    return 1;
  }
  // Setup underlying TriangleMatrix
  Setup( Ntemp );
  // Read number of elements (just another check really).
  // TODO: read long int?
  infile.Read( &Ntemp, sizeof(int) );
  if ( Ntemp != (int)nelements_ ) {
    mprinterr("Internal Error: ClusterMatrix setup unsuccessful!\n");
    return 1;
  }
  // Read elements
  infile.Read( elements_, Ntemp*sizeof(float) );
  // TODO: Read ignore?
  SetupIgnore();
  return 0;
}

// ClusterMatrix::SetupIgnore()
int ClusterMatrix::SetupIgnore() {
  ignore_.assign( nrows_, false );
  return 0;
}

// ClusterMatrix::FindMin()
/** Find the minimum; set corresponding row and column. */
double ClusterMatrix::FindMin(int& iOut, int& jOut) const {
  iOut = -1;
  jOut = -1;
  int iVal = 0;
  int jVal = 1;
  float min = FLT_MAX;
  for (size_t idx = 0; idx < nelements_; ++idx) {
    if (ignore_[iVal] || ignore_[jVal]) {
      // If we dont care about this row/col, just increment
      jVal++;
      if (jVal == nrows_) {
        iVal++;
        jVal = iVal + 1;
      }
    } else {
      // Otherwise search for minimum
      if ( elements_[idx] < min ) {
        min = elements_[idx];
        iOut = iVal;
        jOut = jVal;
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

// ClusterMatrix::PrintElements()
void ClusterMatrix::PrintElements() const {
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
