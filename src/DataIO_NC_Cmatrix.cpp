#include "DataIO_NC_Cmatrix.h"
#include "DataSet_Cmatrix_MEM.h"

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

  return 0;
}
