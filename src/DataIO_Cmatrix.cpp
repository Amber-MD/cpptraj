#include "DataIO_Cmatrix.h"
#include "CpptrajStdio.h" 

// NOTES:
//   Version 1: Add write of ignore array when reduced. Write nrows and
//              and nelements as 8 byte integers.
//   Version 2: Instead of nrows and nelements, write original nrows
//              and actual nrows to easily determine if this is a reduced
//              matrix. Also write sieve value.
//   Version 2 Update: Read/write sieve value as signed, negative
//                     value is random sieve. Variable is same #
//                     of bytes so should be backwards-compatible.
const unsigned char DataIO_Cmatrix::Magic_[4] = {'C', 'T', 'M', 2};

// CONSTRUCTOR
DataIO_Cmatrix::DataIO_Cmatrix()
{
  SetValid( DataSet::CMATRIX );
}

bool DataIO_Cmatrix::ID_DataFormat(CpptrajFile& infile) {
  unsigned char magic[4];
  if (infile.OpenFile()) return false;
  infile.Read( magic, 4 );
  infile.CloseFile();
  return (magic[0]==Magic_[0] && magic[1]==Magic_[1] && magic[2]==Magic_[2]);
}

// -----------------------------------------------------------------------------
// DataIO_Cmatrix::ReadHelp()
void DataIO_Cmatrix::ReadHelp() {

}

// DataIO_Cmatrix::processReadArgs()
int DataIO_Cmatrix::processReadArgs(ArgList& argIn) {

  return 0;
}

// DataIO_Cmatrix::ReadData() 
int DataIO_Cmatrix::ReadData(FileName const& fname, 
                             DataSetList& dsl, std::string const& dsname)
{
  // Allocate data set
  MetaData md( dsname, MetaData::M_MATRIX );
  DataSet* ds = dsl.AddSet(DataSet::CMATRIX, md, "Cmatrix");
  if (ds == 0) return 1;
  DataSet_Cmatrix& Mat = static_cast<DataSet_Cmatrix&>( *ds );
  return ReadCmatrix(fname, Mat);
}

// DataIO_Cmatrix::ReadCmatrix()
int DataIO_Cmatrix::ReadCmatrix(FileName const& fname, DataSet_Cmatrix& Mat) {
  unsigned char magic[4];
  CpptrajFile infile;
  uint_8 ROWS, ELTS;
  sint_8 SIEVE;
  int sieve = 1;
  size_t actual_nrows = 0;
  // Open file for reading
  if (infile.OpenRead(fname)) {
    mprinterr("Error: Could not open '%s' for read.\n", fname.full());
    return 1;
  }
  // SANITY CHECK: Read and check magic byte
  infile.Read( magic, 4 );
  if ( magic[0]!=Magic_[0] || magic[1]!=Magic_[1] || magic[2]!=Magic_[2] ) {
    mprinterr("Error: File '%s' is not a Cpptraj Cluster Matrix file.\n", fname.full());
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
    infile.Read( &SIEVE, sizeof(sint_8) ); // V2: Sieve
    sieve = (int)SIEVE;
  } else {
    mprinterr("Error: ClusterMatrix version %u is not recognized.\n", (unsigned int)magic[3]);
    return 1;
  }
  // If number of rows is not what was expected, abort TODO reimplement somewhere else
/*
  if (ROWS != (uint_8)sizeIn) {
    mprinterr("Error: ClusterMatrix file %s has %lu rows, expected %i.\n",
              fname.full(), ROWS, sizeIn);
    return 1;
  }
*/
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
    sieve = 1;
  }
  // Setup underlying TriangleMatrix for actual # of rows
  if ( Mat.AllocateTriangle( actual_nrows ) ) return 1;
  // Read in matrix elements
  infile.Read( Mat.Ptr(), Mat.Size()*sizeof(float) );
  // If sieved, read in the ignore array
  std::vector<char> ignore_in;
  if (sieve != 1) {
    mprintf("Warning: ClusterMatrix %s contains sieved data.\n", fname.full());
    ignore_in.resize( ROWS ); // Original nrows
    infile.Read( &ignore_in[0], ROWS*sizeof(char) );
  }
  // Setup ignore array and sieve; if not sieving all elements set to false.
  if (Mat.SetupIgnore(ROWS, ignore_in, sieve)) return 1;

  return 0;
}

// -----------------------------------------------------------------------------
// DataIO_Cmatrix::WriteHelp()
void DataIO_Cmatrix::WriteHelp() {

}

// DataIO_Cmatrix::processWriteArgs()
int DataIO_Cmatrix::processWriteArgs(ArgList &argIn) {

  return 0;
}

// DataIO_Cmatrix::WriteData()
int DataIO_Cmatrix::WriteData(FileName const& fname, DataSetList const& SetList)
{
  if (SetList.empty()) return 1;
  if (SetList.size() > 1)
    mprintf("Warning: Multiple sets not yet supported for cluster matrix write.\n");
  DataSet_Cmatrix const& Mat = static_cast<DataSet_Cmatrix const&>( *(*(SetList.begin())) );
  return WriteCmatrix( fname, Mat );
}

// DataIO_Cmatrix::WriteCmatrix()
int DataIO_Cmatrix::WriteCmatrix(FileName const& fname, DataSet_Cmatrix const& Mat) {
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
  // Write magic byte
  outfile.Write( Magic_, 4 );
  // Write original nrows (size of ignore)
  ntemp = (uint_8)Mat.Nframes();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write actual nrows
  ntemp = (uint_8)Mat.Nrows();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write out sieve value
  sint_8 stemp = (sint_8)Mat.SieveValue();
  outfile.Write( &stemp, sizeof(sint_8) );
  // Write matrix elements
  outfile.Write( Mat.Ptr(), Mat.Size()*sizeof(float) );
  // If this is a reduced matrix, write the ignore array as chars.
  if (Mat.SieveType() != ClusterSieve::NONE) {
    std::vector<char> ignore_out( Mat.Nframes() );
    for (unsigned int idx = 0; idx != Mat.Nframes(); idx++) 
      if (Mat.IgnoringRow(idx))
        ignore_out[idx] = 'T';
      else
        ignore_out[idx] = 'F';
    outfile.Write( &ignore_out[0], Mat.Nframes()*sizeof(char) );
  }
  return 0;
}
