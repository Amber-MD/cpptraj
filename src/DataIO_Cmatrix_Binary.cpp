#include "DataIO_Cmatrix_Binary.h"
#include "CpptrajStdio.h"
#include "Cluster/Cmatrix_Binary.h"

using namespace Cpptraj;
using namespace Cluster; 

// CONSTRUCTOR
DataIO_Cmatrix_Binary::DataIO_Cmatrix_Binary()
{
  SetValid( DataSet::PMATRIX_MEM );
}

bool DataIO_Cmatrix_Binary::ID_DataFormat(CpptrajFile& infile) {
  return Binary_Cmatrix::ID_Cmatrix( infile );
}

// -----------------------------------------------------------------------------
// DataIO_Cmatrix_Binary::ReadHelp()
void DataIO_Cmatrix_Binary::ReadHelp() {

}

// DataIO_Cmatrix_Binary::processReadArgs()
int DataIO_Cmatrix_Binary::processReadArgs(ArgList& argIn) {

  return 0;
}

// DataIO_Cmatrix_Binary::ReadData() 
int DataIO_Cmatrix_Binary::ReadData(FileName const& fname, 
                             DataSetList& dsl, std::string const& dsname)
{
  // Allocate data set
  MetaData md( dsname, MetaData::M_MATRIX );
  DataSet* ds = dsl.AddSet(DataSet::PMATRIX_MEM, md, "Cmatrix");
  if (ds == 0) return 1;
  DataSet_PairwiseCache_MEM& Mat = static_cast<DataSet_PairwiseCache_MEM&>( *ds );
  return ReadCmatrix(fname, Mat);
}

// DataIO_Cmatrix_Binary::ReadCmatrix()
int DataIO_Cmatrix_Binary::ReadCmatrix(FileName const& fname, DataSet_PairwiseCache_MEM& Mat) {
  Binary_Cmatrix infile;

  if (infile.OpenCmatrixRead( fname )) return 1;

  // Setup underlying TriangleMatrix for actual # of rows
  if ( Mat.Allocate( DataSet::SizeArray(1, infile.ActualNrows()) ) ) return 1;
  // Allocate sieve status array if needed
  char* sieveStatus = 0;
  std::vector<char> ssArray; // TODO just pass in char array?
  if (infile.Sieve() != 1) {
    ssArray.resize( infile.Ntotal() );
    sieveStatus = &ssArray[0];
  }
  // Read in matrix elements
  if (infile.GetCmatrix( Mat.Ptr(), sieveStatus )) return 1;

  // Set sieve status.
//  if (Mat.SetSieveFromArray(sieveStatus, sieve)) return 1;

  return 0;
}

// -----------------------------------------------------------------------------
/*
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
  DataSet_Cmatrix_MEM const& Mat = static_cast<DataSet_Cmatrix_MEM const&>( *(*(SetList.begin())) );
  return WriteCmatrix( fname, Mat );
}

// DataIO_Cmatrix::WriteCmatrix()
int DataIO_Cmatrix::WriteCmatrix(FileName const& fname, DataSet_Cmatrix_MEM const& Mat) {
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
  // Write original number of frames.
  ntemp = (uint_8)Mat.OriginalNframes();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write actual nrows
  ntemp = (uint_8)Mat.Nrows();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write out sieve value
  sint_8 stemp = (sint_8)Mat.SieveValue();
  outfile.Write( &stemp, sizeof(sint_8) );
  // Write matrix elements
  outfile.Write( Mat.Ptr(), Mat.Size()*sizeof(float) );
  // If this is a reduced matrix, write whether each frame was sieved (T) or not (F). 
  if (Mat.SieveType() != ClusterSieve::NONE) {
    std::vector<char> sieveStatus( Mat.OriginalNframes() );
    for (int idx = 0; idx != Mat.OriginalNframes(); idx++) 
      if (Mat.FrameWasSieved(idx))
        sieveStatus[idx] = 'T';
      else
        sieveStatus[idx] = 'F';
    outfile.Write( &sieveStatus[0], Mat.OriginalNframes()*sizeof(char) );
  }
  return 0;
}
*/
