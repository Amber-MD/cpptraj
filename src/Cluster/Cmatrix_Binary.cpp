#include "Cmatrix_Binary.h"
#include "../CpptrajStdio.h"

Cpptraj::Cluster::Cmatrix_Binary::Cmatrix_Binary() :
  sieve_(0),
  actual_nrows_(0),
  ntotal_(0),
  headerOffset_(0) 
{}

// NOTES:
//   Version 1: Add write of ignore array when reduced. Write nrows and
//              and nelements as 8 byte integers.
//   Version 2: Instead of nrows and nelements, write original nrows
//              and actual nrows to easily determine if this is a reduced
//              matrix. Also write sieve value.
//   Version 2 Update: Read/write sieve value as signed, negative
//                     value is random sieve. Variable is same #
//                     of bytes so should be backwards-compatible.
const unsigned char Cpptraj::Cluster::Cmatrix_Binary::Magic_[4] = {'C', 'T', 'M', 2};

/** File must be set up read. */
bool Cpptraj::Cluster::Cmatrix_Binary::ID_Cmatrix(CpptrajFile& infile)
{
  unsigned char magic[4];
  
  if (infile.OpenFile()) return false;
  infile.Read( magic, 4 );
  infile.CloseFile();
  return (magic[0]==Magic_[0] && magic[1]==Magic_[1] && magic[2]==Magic_[2]);
}

int Cpptraj::Cluster::Cmatrix_Binary::OpenCmatrixRead(FileName const& fname)
{
  unsigned char magic[4];
  uint_8 ROWS, ELTS;
  sint_8 SIEVE;
  sieve_ = 1;
  actual_nrows_ = 0;

  // Open file for reading
  if (file_.OpenRead(fname)) {
    mprinterr("Error: Could not open '%s' for read.\n", fname.full());
    return 1;
  }
  // SANITY CHECK: Read and check magic byte
  file_.Read( magic, 4 );
  if ( magic[0]!=Magic_[0] || magic[1]!=Magic_[1] || magic[2]!=Magic_[2] ) {
    mprinterr("Error: File '%s' is not a Cpptraj Cluster Matrix file.\n", fname.full());
    return 1;
  }
  // Check version, read in nrows and nelements.
  if (magic[3] == 0) {
    int Ntemp = 0;
    file_.Read( &Ntemp, sizeof(int) );
    ROWS = (uint_8)Ntemp;
    actual_nrows_ = (size_t)ROWS;
    file_.Read( &Ntemp, sizeof(int) );
    ELTS = (uint_8)Ntemp;
  } else if (magic[3] == 1) {
    file_.Read( &ROWS, sizeof(uint_8) );
    actual_nrows_ = (size_t)ROWS;
    file_.Read( &ELTS, sizeof(uint_8) );
  } else if (magic[3] == 2) {
    file_.Read( &ROWS,  sizeof(uint_8) ); // V2: Original Nrows
    file_.Read( &ELTS,  sizeof(uint_8) ); // V2: Actual Nrows
    actual_nrows_ = (size_t)ELTS;
    file_.Read( &SIEVE, sizeof(sint_8) ); // V2: Sieve
    sieve_ = (int)SIEVE;
  } else {
    mprinterr("Error: ClusterMatrix version %u is not recognized.\n", (unsigned int)magic[3]);
    return 1;
  }
  // Save current header position
  headerOffset_ = file_.Tell();
  // Checks
  if (magic[3] == 0 || magic[3] == 1) {
    // Version 0/1: Actual # of rows is not known yet. Check that the # elements
    // in the file match the original # elements (i.e. matrix is not sieved).
    // If it is sieved this is not supported.
    uint_8 original_nelements = ( ROWS * (ROWS - 1UL) ) / 2UL;
    if ( original_nelements != ELTS ) {
      mprinterr("Error: Sieved data in ClusterMatrix file %s (version %u) not supported.\n",
                fname.full(), (unsigned int)magic[3]);
      return 1;
    }
    sieve_ = 1;
  }
  ntotal_ = (size_t)ROWS;
  return 0;
}

static long int calcTriIndex(size_t nX, size_t xIn, size_t yIn) {
      size_t i, j;
      if (yIn > xIn) {
        i = xIn;
        j = yIn;
      } else if (xIn > yIn) {
        i = yIn;
        j = xIn;
      } else // iIn == jIn, triangle matrix diagonal is indicated by -1 
        return -1L;
      size_t i1 = i + 1UL;
      return (long int)(( (nX * i) - ((i1 * i) / 2UL) ) + j - i1);
}

double Cpptraj::Cluster::Cmatrix_Binary::GetCmatrixElement(unsigned int col, unsigned int row)
{
  // Determine absolute index
  long int idx = calcTriIndex(actual_nrows_, col, row);
  file_.Seek( headerOffset_ + idx * sizeof(float) );
  float fvar;
  file_.Read( &fvar, sizeof(float) );
  return (double)fvar;
}

double Cpptraj::Cluster::Cmatrix_Binary::GetCmatrixElement(unsigned int idx)
{
  file_.Seek( headerOffset_ + idx * sizeof(float) );
  float fvar;
  file_.Read( &fvar, sizeof(float) );
  return (double)fvar;
}

static inline size_t Nelements(size_t nrows) {
  return ( nrows * (nrows - 1UL) ) / 2UL;
}

int Cpptraj::Cluster::Cmatrix_Binary::GetCmatrix(float* ptr, char* frameIsPresent) {
  file_.Seek( headerOffset_  );
  size_t nelements =  Nelements( actual_nrows_ );
  if (file_.Read( ptr, nelements * sizeof(float) ) < 1) return 1;
  // Read sieve if needed
  if (frameIsPresent != 0) {
    if (file_.Read( frameIsPresent, ntotal_ * sizeof(char) ) < 1) return 1;
  }
  return 0;
}

int Cpptraj::Cluster::Cmatrix_Binary::WriteCmatrix(FileName const& fname,
                                                   const float* ptr,
                                                   Cframes const& frameToIdx,
                                                   size_t Nrows,
                                                   int sieveValue)
{
  CpptrajFile outfile;
  uint_8 ntemp;
  // No stdout write allowed.
  if (fname.empty()) {
    mprinterr("Internal Error: DataIO_Cmatrix::WriteData() called with no filename.\n");
    return 1;
  }
  if (outfile.OpenWrite(fname)) {
    mprinterr("Error: Could not open %s for write.\n", fname.full());
    return 1;
  }
  size_t OriginalNframes = frameToIdx.size();
  // Write magic byte
  outfile.Write( Magic_, 4 );
  // Write original number of frames.
  ntemp = (uint_8)OriginalNframes;
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write actual nrows
  ntemp = (uint_8)Nrows;
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write out sieve value
  sint_8 stemp = (sint_8)sieveValue;
  outfile.Write( &stemp, sizeof(sint_8) );
  // Write matrix elements
  size_t nelements =  Nelements( Nrows );
  outfile.Write( ptr, nelements*sizeof(float) );
  // If this is a reduced matrix, write whether each frame was sieved (T) or not (F). 
  if (sieveValue != 1) {
    //DataSet_PairwiseMatrix::StatusArray frameIsPresent;
    std::vector<char> frameIsPresent;
    frameIsPresent.reserve( OriginalNframes );
    for (Cframes::const_iterator it = frameToIdx.begin();
                                 it != frameToIdx.end(); ++it)
      if (*it == -1)
        frameIsPresent.push_back( 'F' );
      else
        frameIsPresent.push_back( 'T' );
    outfile.Write( &frameIsPresent[0], OriginalNframes*sizeof(char) );
  }
  return 0;
}
