#include <cfloat> // FLT_MAX
#include "ClusterMatrix.h"
#include "CpptrajStdio.h"

// NOTES:
//   Version 1: Add write of ignore array when reduced. Write nrows and
//              and nelements as 8 byte integers.
const unsigned char ClusterMatrix::Magic_[4] = {'C', 'T', 'M', 1};

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
  uint_8 ntemp;
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
  // Determine if any rows are ignored. Will write reduced matrix if so.
  int nignored = 0;
  for (std::vector<bool>::const_iterator ig = ignore_.begin(); ig != ignore_.end(); ++ig) {
    if (*ig) ++nignored; 
  }
  if (nignored > 0) {
    mprintf("\tSaving reduced matrix with %i rows/cols\n", Nrows() - nignored);
    //mprintf("Ignoring %i out of %i rows/cols\n", nignored, Nrows()); // DEBUG
    ClusterMatrix reduced(Nrows() - nignored);
    for (int row = 0; row < Nrows() - 1; ++row) {
      if (!ignore_[row]) {
        for (int col = row + 1; col < Nrows(); ++col) {
          if (!ignore_[col]) 
            reduced.AddElement( this->GetElementF(row, col) );
        }
      }
    }
    // DEBUG - Print reduced matrix
    //mprintf("Reduced matrix:\n");
    //reduced.PrintElements();
    // Write original nrows
    ntemp = (uint_8)nrows_;
    outfile.Write( &ntemp, sizeof(uint_8) );
    // Write reduced # of elements
    ntemp = (uint_8)reduced.nelements_;
    outfile.Write( &ntemp, sizeof(uint_8) );
    // Write reduced matrix elements
    outfile.Write( reduced.elements_, reduced.nelements_*sizeof(float) );
    // Write ignore array - convert to char array first
    char* ignore_out = new char[ ignore_.size() ];
    int idx = 0;
    for (std::vector<bool>::const_iterator ig = ignore_.begin(); ig != ignore_.end(); ++ig) 
      if (*ig)
        ignore_out[idx++] = 'T';
      else
        ignore_out[idx++] = 'F';
    outfile.Write( ignore_out, ignore_.size()*sizeof(char) );
    delete[] ignore_out;
  } else {
    // Write nrows
    ntemp = (uint_8)nrows_;
    outfile.Write( &ntemp, sizeof(uint_8) );
    // Write number of elements.
    ntemp = (uint_8)nelements_;
    outfile.Write( &ntemp, sizeof(uint_8) );
    // Write matrix elements
    outfile.Write( elements_, nelements_*sizeof(float) );
  }
  return 0;
}

// ClusterMatrix::LoadFile()
int ClusterMatrix::LoadFile(std::string const& filename, int sizeIn) {
  unsigned char magic[4];
  CpptrajFile infile;
  uint_8 ROWS;
  uint_8 ELTS;
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
  // Check version, read in nrows and nelements.
  if (magic[3] == 0) {
    int Ntemp = 0;
    infile.Read( &Ntemp, sizeof(int) );
    ROWS = (uint_8)Ntemp;
    infile.Read( &Ntemp, sizeof(int) );
    ELTS = (uint_8)Ntemp;
  } else if (magic[3] == 1) {
    infile.Read( &ROWS, sizeof(uint_8) );
    infile.Read( &ELTS, sizeof(uint_8) );
  } else {
    mprinterr("Error: ClusterMatrix version %u is not recognized.\n", (unsigned int)magic[3]);
    return 1;
  }
  // If number of rows is not what was expected, abort
  if (ROWS != (uint_8)sizeIn) {
    mprinterr("Error: ClusterMatrix file %s has %lu rows, expected %i.\n",
              filename.c_str(), ROWS, sizeIn);
    return 1;
  }
  // Setup underlying TriangleMatrix
  Setup( ROWS );
  SetupIgnore();
  // Read in matrix
  if ( nelements_ != (size_t)ELTS ) {
    if ( magic[3] == 0 ) {
      mprinterr("Error: ClusterMatrix %s is version 0, does not support sieved data.\n",
                filename.c_str());
      return 1;
    }
    mprintf("Warning: ClusterMatrix %s contains sieved data.\n", filename.c_str());
    mprintf("\tReading in reduced matrix with %lu elements.\n", ELTS);
    // Reduced matrix. Need to read elements and ignore array.
    float* reduced = new float[ ELTS ];
    infile.Read( reduced, ELTS*sizeof(float) );
    float* RF = reduced;
    char* ignore_in = new char[ ROWS ];
    infile.Read( ignore_in, ROWS*sizeof(char) );
    for (int row = 0; row < Nrows() - 1; ++row) {
      if (ignore_in[row] == 'F') {
        for (int col = row + 1; col < Nrows(); ++col) {
          if (ignore_in[col] == 'F') 
            SetElementF( row, col, *(RF++) );
        }
      } else {
        ignore_[row] = true;
      }
    }
    delete[] ignore_in;
    delete[] reduced;
    //PrintElements(); // DEBUG
  } else {
    // Complete matrix. Read elements only.
    infile.Read( elements_, ELTS*sizeof(float) );
  }
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
  size_t iVal = 0L;
  size_t jVal = 1L;
  float min = FLT_MAX;
  for (size_t idx = 0L; idx < nelements_; ++idx) {
    if (ignore_[iVal] || ignore_[jVal]) {
      // If we dont care about this row/col, just increment
      jVal++;
      if (jVal == nrows_) {
        iVal++;
        jVal = iVal + 1L;
      }
    } else {
      // Otherwise search for minimum
      if ( elements_[idx] < min ) {
        min = elements_[idx];
        iOut = (int)iVal;
        jOut = (int)jVal;
      }
      // Increment indices
      jVal++;
      if (jVal == nrows_) {
        iVal++;
        jVal = iVal + 1L;
      }
    }
  }
  return (double)min;
}

// ClusterMatrix::PrintElements()
void ClusterMatrix::PrintElements() const {
  size_t iVal = 0;
  size_t jVal = 1;
  for (size_t idx = 0L; idx < nelements_; idx++) {
    if (!ignore_[iVal] && !ignore_[jVal])
      mprintf("\t%u %u %8.3f\n",iVal,jVal,elements_[idx]);
    // Increment indices
    jVal++;
    if (jVal == nrows_) {
      ++iVal;
      jVal = iVal + 1L;
    }
  }
} 
