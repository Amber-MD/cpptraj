#include "Matrix_2D.h"

// CONSTRUCTOR
Matrix_2D::Matrix_2D() : DataSet(MATRIX2D, 12, 4, 2) {} 

void Matrix_2D::Write2D(CpptrajFile& outfile, int xIn, int yIn) const {
  size_t x = (size_t)xIn;
  size_t y = (size_t)yIn;
  if ( x < 0 || y < 0 || x >= mat_.Ncols() || y >= mat_.Nrows() )
    outfile.Printf(data_format_, 0.0);
  else 
    outfile.Printf(data_format_, mat_.element(x,y));
}
