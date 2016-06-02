#include "DataSet_Cmatrix.h"
#include "CpptrajStdio.h"

void DataSet_Cmatrix::PrintElements() const {
  // NOTE: Matrix is square upper triangle, Nrows == Ncols
  for (unsigned int row = 0; row != Nrows(); row++)
    for (unsigned int col = row + 1; col != Nrows(); col++)
      mprintf("\t%u %u %f\n", row+1, col+1, GetFdist(col, row));
}
