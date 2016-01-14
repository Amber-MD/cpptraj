#include "DataSet_MatrixFlt.h"
void DataSet_MatrixFlt::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = (size_t)pIn[0];
  size_t y = (size_t)pIn[1];
  if ( x >= mat_.Ncols() || y >= mat_.Nrows() )
    outfile.Printf(format_.fmt(), 0.0);
  else 
    outfile.Printf(format_.fmt(), mat_.element(x,y));
}

double* DataSet_MatrixFlt::MatrixArray() const {
  double* matOut = new double[ mat_.size() ];
  for (size_t i = 0; i < mat_.size(); ++i)
    matOut[i] = (double)mat_[i];
  return matOut;
}
#ifdef MPI
int DataSet_MatrixFlt::Sync(size_t total, std::vector<int> const& rank_frames,
                            Parallel::Comm const& commIn)
{
  if (commIn.Master()) {
    std::vector<float> buf( mat_.size() );
    commIn.Reduce( &(buf[0]), &(mat_[0]),  mat_.size(),  MPI_FLOAT, MPI_SUM );
    std::copy( buf.begin(), buf.end(), mat_.begin() );
  } else
    commIn.Reduce( 0,         &(mat_[0]),  mat_.size(),  MPI_FLOAT, MPI_SUM );
  return 0;
}
#endif
