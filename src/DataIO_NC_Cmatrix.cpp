#include "DataIO_NC_Cmatrix.h"
#include "DataSet_Cmatrix_MEM.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataIO_NC_Cmatrix::DataIO_NC_Cmatrix() {
  SetValid( DataSet::CMATRIX );
}

// DataIO_NC_Cmatrix::ReadData()
int DataIO_NC_Cmatrix::ReadData(FileName const& fname,
                         DataSetList& dsl, std::string const& dsname)
{
  int sieve;
  if (file_.OpenCmatrixRead(fname, sieve)) return 1;

  MetaData md( dsname );
  md.SetFileName( fname );
  DataSet* ds = dsl.AddSet( DataSet::CMATRIX, md );
  if (ds == 0) return 1;
  DataSet_Cmatrix_MEM& cmatrix = static_cast<DataSet_Cmatrix_MEM&>( *ds );
  // Allocate matrix for actual number of rows
  if (cmatrix.Allocate( DataSet::SizeArray(1, file_.MatrixRows()) )) return 1;
  // Get the sieve status array.
  if (cmatrix.SetSieveFromArray( file_.GetSieveStatus(), sieve )) return 1;
  // Get the matrix
  if (file_.GetCmatrix( cmatrix.Ptr() )) return 1;
  file_.CloseCmatrix();

  return 0;
}

int DataIO_NC_Cmatrix::WriteData(FileName const& fname, DataSetList const& SetList)
{
  if (SetList.empty()) return 1;
  if (SetList.size() > 1)
    mprintf("Warning: Multiple sets not yet supported for cluster matrix write.\n");
  DataSet_Cmatrix_MEM const& Mat = static_cast<DataSet_Cmatrix_MEM const&>( *(*(SetList.begin())) );
  // Create the file
  if (file_.CreateCmatrix( fname, Mat.OriginalNframes(), Mat.Nrows(), Mat.SieveValue(),
                           Mat.MetricDescription() ))
    return 1;
  // Write the matrix
  if (file_.WriteCmatrix( Mat.Ptr() )) return 1;
  // Write sieved frames
  if (Mat.SieveType() != ClusterSieve::NONE) {
    if (file_.WriteFramesArray( Mat.FramesToCluster() )) return 1;
  }

  file_.CloseCmatrix();
  return 0;
}
