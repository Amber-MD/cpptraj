#include <cfloat> // FLT_MAX
#include "ClusterMatrix.h"
#include "CpptrajStdio.h"

// NOTES:
//   Version 1: Add write of ignore array when reduced. Write nrows and
//              and nelements as 8 byte integers.
//   Version 2: Instead of nrows and nelements, write original nrows
//              and actual nrows to easily determine if this is a reduced
//              matrix. Also write sieve value.
const unsigned char ClusterMatrix::Magic_[4] = {'C', 'T', 'M', 2};

// CONSTRUCTOR
/** Intended for use with cluster pairwise distance calculations
  * where frames may be sieved. The underlying TriangleMatrix will
  * only be set up to hold the actual number of frames based on
  * the sieve value. The Ignore array will be set up based on
  * the original number of frames.
  */
ClusterMatrix::ClusterMatrix(size_t sizeIn, size_t sieveIn) : sieve_(sieveIn)
{
  if (sieve_ < 1) sieve_ = 1;
  // Sieved distances should be ignored
  if (sieve_ > 1) {
    // Determine size of underlying TriangleMatrix
    size_t actual_nrows = sizeIn / sieve_;
    if ( (sizeIn % sieve_) > 0 )
      ++actual_nrows;
    TriangleMatrix::Setup( actual_nrows );
    // Set up the ignore array to ignore sieved frames
    ignore_.assign(sizeIn, true);
    for (size_t frame = 0; frame < sizeIn; frame += sieve_) 
      ignore_[frame] = false;
  } else {
    TriangleMatrix::Setup( sizeIn );
    ignore_.assign(sizeIn, false);
  }
}

// COPY CONSTRUCTOR
ClusterMatrix::ClusterMatrix(const ClusterMatrix& rhs) :
  TriangleMatrix(rhs),
  ignore_(rhs.ignore_),
  sieve_(rhs.sieve_)
{}

// ASSIGNMENT
ClusterMatrix& ClusterMatrix::operator=(const ClusterMatrix& rhs) {
  if (this == &rhs) return *this;
  TriangleMatrix::operator=(rhs);
  ignore_ = rhs.ignore_;
  sieve_ = rhs.sieve_;
  return *this;
}

// ClusterMatrix::SaveFile()
/** Save the matrix to a binary file. Format is:
  *     Full: [4*char]     [uint_8] [uint_8]    [uint_8] [nelements*float]
  *  Reduced: [4*char]     [uint_8] [uint_8]    [uint_8] [nelements*float] [nrows*char]
  *     Vars: ['C''T''M'X] [nrows]  [nelements] [sieve]  [elements]        [ignore]
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
  // Write original nrows (size of ignore)
  ntemp = (uint_8)ignore_.size();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write actual nrows
  ntemp = (uint_8)Nrows();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write out sieve value
  ntemp = (uint_8)sieve_;
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write matrix elements
  outfile.Write( elements_, Nelements()*sizeof(float) );
  // If this is a reduced matrix, write the ignore array as chars.
  if (sieve_ > 1) {
    char* ignore_out = new char[ ignore_.size() ];
    int idx = 0;
    for (std::vector<bool>::const_iterator ig = ignore_.begin(); ig != ignore_.end(); ++ig) 
      if (*ig)
        ignore_out[idx++] = 'T';
      else
        ignore_out[idx++] = 'F';
    outfile.Write( ignore_out, ignore_.size()*sizeof(char) );
    delete[] ignore_out;
  }
  return 0;
}

// ClusterMatrix::LoadFile()
int ClusterMatrix::LoadFile(std::string const& filename, int sizeIn) {
  unsigned char magic[4];
  CpptrajFile infile;
  uint_8 ROWS, ELTS, SIEVE;
  size_t actual_nrows = 0;
  // Open file for reading
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
    actual_nrows = (size_t)ROWS;
    infile.Read( &Ntemp, sizeof(int) );
    ELTS = (uint_8)Ntemp;
  } else if (magic[3] == 1) {
    infile.Read( &ROWS, sizeof(uint_8) );
    actual_nrows = (size_t)ROWS;
    infile.Read( &ELTS, sizeof(uint_8) );
  } else if (magic[3] == 2) {
    infile.Read( &ROWS,  sizeof(uint_8) ); // V2: Original Nrows
    infile.Read( &ELTS,  sizeof(uint_8) ); // V2: Actual Nrows
    actual_nrows = (size_t)ELTS;
    infile.Read( &SIEVE, sizeof(uint_8) ); // V2: Sieve
    sieve_ = (size_t)SIEVE;
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
  if (magic[3] == 0 || magic[3] == 1) {
    // Version 0/1: Actual # of rows is not known yet. Check that the # elements
    // in the file match the original # elements (i.e. matrix is not sieved).
    // If it is sieved this is not supported.
    uint_8 original_nelements = ( ROWS * (ROWS - 1UL) ) / 2UL;
    if ( original_nelements != ELTS ) {
      mprinterr("Error: Sieved data in ClusterMatrix file %s (version %u) not supported.\n",
                filename.c_str(), (unsigned int)magic[3]);
      return 1;
    }
    sieve_ = 1;
  }
  // Setup underlying TriangleMatrix for actual # of rows
  if (TriangleMatrix::Setup( actual_nrows )) return 1;
  // Set all ignore elements for original # rows to false.
  ignore_.assign(ROWS, false);
  // Read in matrix elements
  infile.Read( elements_, Nelements()*sizeof(float) );
  // If sieved, read in the ignore array
  if (sieve_ > 1) {
    mprintf("Warning: ClusterMatrix %s contains sieved data.\n", filename.c_str());
    char* ignore_in = new char[ ROWS ]; // Original nrows
    infile.Read( ignore_in, ROWS*sizeof(char) );
    for (uint_8 row = 0; row < ROWS; ++row)
      if (ignore_in[row] == 'T')
        ignore_[row] = true;
    delete[] ignore_in;
  }
  mprintf("\tLoaded %s: %u original rows, %u actual rows, %u elements, sieve=%u\n",
          filename.c_str(), ROWS, Nrows(), Nelements(), sieve_);
  return 0;
}

// ClusterMatrix::SetupMatrix()
int ClusterMatrix::SetupMatrix(size_t sizeIn) {
  if (TriangleMatrix::Setup(sizeIn)) return 1;
  ignore_.assign( sizeIn, false );
  sieve_ = 1;
  return 0;
}

// ClusterMatrix::FindMin()
/** Find the minimum; set corresponding row and column. */
double ClusterMatrix::FindMin(int& iOut, int& jOut) const {
  iOut = -1;
  jOut = -1;
  size_t iVal = 0L;
  size_t jVal = sieve_;
  float min = FLT_MAX;
  for (size_t idx = 0L; idx < Nelements(); ++idx) {
    if (ignore_[iVal] || ignore_[jVal]) {
      // If we dont care about this row/col, just increment
      jVal += sieve_;
      if (jVal >= ignore_.size()) {
        iVal += sieve_;
        jVal = iVal + sieve_;
      }
    } else {
      // Otherwise search for minimum
      if ( elements_[idx] < min ) {
        min = elements_[idx];
        iOut = (int)iVal;
        jOut = (int)jVal;
      }
      // Increment indices
      jVal += sieve_;
      if (jVal >= ignore_.size()) {
        iVal += sieve_;
        jVal = iVal + sieve_;
      }
    }
  }
  return (double)min;
}

// ClusterMatrix::PrintElements()
void ClusterMatrix::PrintElements() const {
  size_t iVal = 0L;
  size_t jVal = sieve_;
  for (size_t idx = 0L; idx < Nelements(); idx++) {
    if (!ignore_[iVal] && !ignore_[jVal])
      mprintf("\t%u %u %8.3f\n",iVal,jVal,elements_[idx]);
    // Increment indices
    jVal += sieve_;
    if (jVal >= ignore_.size()) {
      iVal += sieve_;
      jVal = iVal + sieve_;
    }
  }
}
