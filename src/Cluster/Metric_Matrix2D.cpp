#include "Metric_Matrix2D.h"

static inline void idxToColRow(int idx, int Ncols, int& col, int& row)
{
  row = idx / Ncols;
  col = idx - (row * Ncols);
}


double Cpptraj::Cluster::Metric_Matrix2D::FrameDist(int f1, int f2)
{
  //int col, row;
  //idxToColRow( f1, matrix_->Ncols(), col, row );
  double val1 = matrix_->GetElement( f1 );
  double val2 = matrix_->GetElement( f2 );
  int col1, row1;
  idxToColRow( f1, matrix_->Ncols(), col1, row1 );
  int col2, row2;
  idxToColRow( f2, matrix_->Ncols(), col2, row2 );

  double dv = val1 - val2;
  double dr = (double)(row1 - row2);
  double dc = (double)(col1 - col2);
  double dist2 = dv*dv + dr*dr + dc*dc;
  return dist2;
}
