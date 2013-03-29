#include "DataSet_MatrixVec3.h"
void DataSet_MatrixVec3::Write2D(CpptrajFile& outfile, int xIn, int yIn) const {
  size_t x = (size_t)xIn;
  size_t y = (size_t)yIn;
  if ( x < 0 || y < 0 || x >= mat_.Ncols() || y >= mat_.Nrows() )
    outfile.Printf(data_format_, 0.0);
  else {
    Vec3 vecOut = mat_.element(x,y); 
    outfile.Printf(data_format_, vecOut[0], vecOut[1], vecOut[2]);
  }
}
