#include "Binary_Cmatrix.h"
#include "../CpptrajStdio.h"

Cpptraj::Cluster::Binary_Cmatrix::Binary_Cmatrix() 
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
const unsigned char Cpptraj::Cluster::Binary_Cmatrix::Magic_[4] = {'C', 'T', 'M', 2};


bool Cpptraj::Cluster::Binary_Cmatrix::ID_Cmatrix(FileName const& fname)
{
  unsigned char magic[4];
  
  CpptrajFile infile;
  if (infile.OpenRead(fname)) return false;
  infile.Read( magic, 4 );
  infile.CloseFile();
  return (magic[0]==Magic_[0] && magic[1]==Magic_[1] && magic[2]==Magic_[2]);
}

int Cpptraj::Cluster::Binary_Cmatrix::OpenCmatrixRead(FileName const& fname)
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
  return 0;
}
