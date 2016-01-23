#include "DataSet_GridFlt.h"
#include "CpptrajStdio.h" // DEBUG

DataSet_GridFlt::DataSet_GridFlt(DataSet_GridFlt const& rhs) : DataSet_3D(rhs), grid_(rhs.grid_) {}

void DataSet_GridFlt::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = pIn[0];
  size_t y = pIn[1];
  size_t z = pIn[2];
  if ( x >= grid_.NX() || y >= grid_.NY() || z >= grid_.NZ() )
    outfile.Printf(format_.fmt(), 0.0);
  else
    outfile.Printf(format_.fmt(), grid_.element(x,y,z));
}

#ifdef MPI
int DataSet_GridFlt::Sync(size_t total, std::vector<int> const& rank_frames,
                          Parallel::Comm const& commIn)
{
  if (commIn.Master()) {
    std::vector<float> buf( grid_.size() );
    commIn.Reduce( &(buf[0]), &(grid_[0]), grid_.size(), MPI_FLOAT, MPI_SUM );
    std::copy( buf.begin(), buf.end(), grid_.begin() );
  } else
    commIn.Reduce( 0,         &(grid_[0]), grid_.size(), MPI_FLOAT, MPI_SUM );
  return 0;
}
#endif
