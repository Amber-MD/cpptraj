#include "DataSet_MatrixDbl.h"
void DataSet_MatrixDbl::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = (size_t)pIn[0];
  size_t y = (size_t)pIn[1];
  if ( x >= mat_.Ncols() || y >= mat_.Nrows() )
    outfile.Printf(format_.fmt(), 0.0);
  else 
    outfile.Printf(format_.fmt(), mat_.element(x,y));
}

double* DataSet_MatrixDbl::MatrixArray() const {
  double* matOut = new double[ mat_.size() ];
  std::copy( mat_.Ptr(), mat_.Ptr() + mat_.size(), matOut );
  return matOut;
}

/** Clear the matrix. */
void DataSet_MatrixDbl::Clear() {
  mat_.clear();
  vect_.clear();
  mass_.clear();
  //kind_;
  snap_ = 0;
}

/** Normalize the matrix. */
void DataSet_MatrixDbl::Normalize(double norm) {
  for (Matrix<double>::iterator it = mat_.begin(); it != mat_.end(); ++it)
    *it *= norm;
}

#ifdef MPI
int DataSet_MatrixDbl::Sync(size_t total, std::vector<int> const& rank_frames,
                            Parallel::Comm const& commIn)
{
  int total_frames = 0;
  int nframes = (int)snap_;
  commIn.ReduceMaster( &total_frames, &nframes, 1, MPI_INT, MPI_SUM );
  if (commIn.Master()) {
    snap_ = (unsigned int)total_frames;
    Darray buf( mat_.size() );
    commIn.ReduceMaster( &(buf[0]), &(mat_[0]),  mat_.size(),  MPI_DOUBLE, MPI_SUM );
    std::copy( buf.begin(), buf.end(), mat_.begin() );
    buf.assign( vect_.size(), 0.0 );
    commIn.ReduceMaster( &(buf[0]), &(vect_[0]), vect_.size(), MPI_DOUBLE, MPI_SUM );
    std::copy( buf.begin(), buf.end(), vect_.begin() );
  } else {
    commIn.ReduceMaster( 0,         &(mat_[0]),  mat_.size(),  MPI_DOUBLE, MPI_SUM );
    commIn.ReduceMaster( 0,         &(vect_[0]), vect_.size(), MPI_DOUBLE, MPI_SUM );
  }
  return 0;
}
#endif
