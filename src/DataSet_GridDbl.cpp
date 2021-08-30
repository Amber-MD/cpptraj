#include "DataSet_GridDbl.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

DataSet_GridDbl::DataSet_GridDbl(DataSet_GridDbl const& rhs) : DataSet_3D(rhs), grid_(rhs.grid_) {}

DataSet_GridDbl& DataSet_GridDbl::operator=(DataSet_GridDbl const& rhs) {
  if (this == &rhs) return *this;
  DataSet_3D::operator=(rhs);
  grid_ = rhs.grid_;
  return *this;
}

void DataSet_GridDbl::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = pIn[0];
  size_t y = pIn[1];
  size_t z = pIn[2];
  if ( x >= grid_.NX() || y >= grid_.NY() || z >= grid_.NZ() )
    outfile.Printf(format_.fmt(), 0.0);
  else
    outfile.Printf(format_.fmt(), grid_.element(x,y,z));
}

#ifdef MPI
int DataSet_GridDbl::SyncGrid(size_t total, std::vector<int> const& rank_frames,
                          Parallel::Comm const& commIn)
{
  if (commIn.Master()) {
    std::vector<double> buf( grid_.size() );
    commIn.ReduceMaster( &(buf[0]), &(grid_[0]), grid_.size(), MPI_FLOAT, MPI_SUM );
    std::copy( buf.begin(), buf.end(), grid_.begin() );
  } else
    commIn.ReduceMaster( 0,         &(grid_[0]), grid_.size(), MPI_FLOAT, MPI_SUM );
  return 0;
}
#endif
