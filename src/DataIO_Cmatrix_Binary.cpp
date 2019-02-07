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
  return Cmatrix_Binary::ID_Cmatrix( infile );
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
  DataSet* ds = dsl.AddSet(DataSet::PMATRIX_MEM, md );
  if (ds == 0) return 1;
  DataSet_PairwiseCache_MEM& Mat = static_cast<DataSet_PairwiseCache_MEM&>( *ds );
  return ReadCmatrix(fname, Mat);
}

// DataIO_Cmatrix_Binary::ReadCmatrix()
int DataIO_Cmatrix_Binary::ReadCmatrix(FileName const& fname, DataSet_PairwiseCache_MEM& Mat) {
  Cmatrix_Binary infile;

  if (infile.OpenCmatrixRead( fname )) return 1;

  // Setup underlying TriangleMatrix for actual # of rows
  if ( Mat.Allocate( DataSet::SizeArray(1, infile.ActualNrows()) ) ) return 1;
  // Allocate sieve status array if needed
  char* sieveStatus = 0;
  DataSet_PairwiseCache::StatusArray ssArray; // TODO just pass in char array instead of ptr?
  if (infile.Sieve() != 1) {
    ssArray.resize( infile.Ntotal() );
    sieveStatus = &ssArray[0];
  }
  // Read in matrix elements
  if (infile.GetCmatrix( Mat.Ptr(), sieveStatus )) return 1;

  // Set sieve status.
  if (Mat.SetupFromStatus(ssArray, infile.Sieve())) return 1;

  return 0;
}

// -----------------------------------------------------------------------------
// DataIO_Cmatrix_Binary::WriteHelp()
void DataIO_Cmatrix_Binary::WriteHelp() {

}

// DataIO_Cmatrix_Binary::processWriteArgs()
int DataIO_Cmatrix_Binary::processWriteArgs(ArgList &argIn) {

  return 0;
}

// DataIO_Cmatrix_Binary::WriteData()
int DataIO_Cmatrix_Binary::WriteData(FileName const& fname, DataSetList const& SetList)
{
  if (SetList.empty()) return 1;
  if (SetList.size() > 1)
    mprintf("Warning: Multiple sets not yet supported for cluster matrix write.\n");
  DataSet_PairwiseCache_MEM const& Mat = 
    static_cast<DataSet_PairwiseCache_MEM const&>( *(*(SetList.begin())) );
  return WriteCmatrix( fname, Mat );
}

// DataIO_Cmatrix_Binary::WriteCmatrix()
int DataIO_Cmatrix_Binary::WriteCmatrix(FileName const& fname,
                                        DataSet_PairwiseCache_MEM const& Mat)
{
  int err = Cmatrix_Binary::WriteCmatrix( fname, Mat.Ptr(), Mat.FrameToIdx(),
                                          Mat.Nrows(), Mat.SieveVal() );
  return err;
}
