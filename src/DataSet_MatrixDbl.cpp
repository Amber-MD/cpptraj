#include "DataSet_MatrixDbl.h"
void DataSet_MatrixDbl::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = (size_t)pIn[0];
  size_t y = (size_t)pIn[1];
  if ( x >= mat_.Ncols() || y >= mat_.Nrows() )
    outfile.Printf(data_format_, 0.0);
  else 
    outfile.Printf(data_format_, mat_.element(x,y));
}

double* DataSet_MatrixDbl::MatrixArray() const {
  double* matOut = new double[ mat_.size() ];
  std::copy( mat_.Ptr(), mat_.Ptr() + mat_.size(), matOut );
  return matOut;
}
