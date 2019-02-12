#include "DataIO_Cmatrix_NC.h"
#include "DataSet_PairwiseCache_MEM.h" // TODO just pairwisecache
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataIO_Cmatrix_NC::DataIO_Cmatrix_NC() {
  SetValid( DataSet::PMATRIX_MEM );
}

// DataIO_Cmatrix_NC::ReadData()
int DataIO_Cmatrix_NC::ReadData(FileName const& fname,
                         DataSetList& dsl, std::string const& dsname)
{
  int sieve;
  if (file_.OpenCmatrixRead(fname, sieve)) return 1;

  MetaData md( dsname );
  md.SetFileName( fname );
  DataSet* ds = dsl.AddSet( DataSet::PMATRIX_MEM, md );
  if (ds == 0) return 1;
  DataSet_PairwiseCache_MEM& cmatrix = static_cast<DataSet_PairwiseCache_MEM&>( *ds );
  // Allocate matrix for actual number of rows
  if (cmatrix.Allocate( DataSet::SizeArray(1, file_.MatrixRows()) )) return 1;
  // Get the sieve status array.
  if (cmatrix.SetupFromStatus( file_.GetSieveStatus(), sieve )) return 1;
  // Get the matrix
  if (file_.GetCmatrix( cmatrix.Ptr() )) return 1;
  file_.CloseCmatrix();

  return 0;
}

int DataIO_Cmatrix_NC::WriteData(FileName const& fname, DataSetList const& SetList)
{
  if (SetList.empty()) return 1;
  if (SetList.size() > 1)
    mprintf("Warning: Multiple sets not yet supported for cluster matrix write.\n");
  DataSet_PairwiseCache_MEM const& Mat =
    static_cast<DataSet_PairwiseCache_MEM const&>( *(*(SetList.begin())) );
  // Create the file
  if (file_.CreateCmatrix( fname, Mat.FrameToIdx().size(), Mat.Nrows(), Mat.SieveVal(),
                           Mat.MetricDescrip() ))
    return 1;
  // Write the matrix
  if (file_.WriteCmatrix( Mat.Ptr() )) return 1;
  // Write present frames
  if (Mat.SieveVal() != 1) {
    if (file_.WriteFramesArray( Mat.PresentFrames() )) return 1;
  }

  file_.CloseCmatrix();
  return 0;
}
