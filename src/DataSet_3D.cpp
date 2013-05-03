#include "DataSet_3D.h"

// DataSet_3D::setOriginAndSpacing()
int DataSet_3D::setOriginAndSpacing(size_t nx, size_t ny, size_t nz,
                                    double ox, double oy, double oz,
                                    double dx, double dy, double dz)
{
  if (nx == 0 || ny == 0 || nz == 0) return 1;
  ox_ = ox;
  oy_ = oy;
  oz_ = oz;
  dx_ = dx;
  dy_ = dy;
  dz_ = dz;
  mx_ = ox_ + ((double)nx * dx_);
  my_ = oy_ + ((double)ny * dy_);
  mz_ = oz_ + ((double)nz * dz_);
  return 0;
}

// Calc_Origin()
/** For even-spaced grids, origin is center - (N/2)*spacing.
  * For odd-spaced grids, origin is center - ((N-1/2)*spacing)+half_spacing
  */
static double Calc_Origin(int N, double D) {
  int odd = N % 2;
  int half = (N - odd) / 2;
  if (odd)
    return -(((double)half * D) + (D * 0.5));
  else
    return -((double)half * D);
}

// DataSet_3D::setCenterAndSpacing()
int DataSet_3D::setCenterAndSpacing(size_t nx, size_t ny, size_t nz,
                                    double cx, double cy, double cz,
                                    double dx, double dy, double dz)
{
  // Calculate origin from center coordinates.
  double ox = cx + Calc_Origin(nx, dx);
  double oy = cy + Calc_Origin(ny, dy);
  double oz = cz + Calc_Origin(nz, dz);
  return setOriginAndSpacing(nx,ny,nz,ox,oy,oz,dx,dy,dz);
}
