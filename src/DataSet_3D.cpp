#include "DataSet_3D.h"
#include "CpptrajStdio.h"

// DataSet_3D::Allocate_N_O_D()
int DataSet_3D::Allocate_N_O_D(size_t nx, size_t ny, size_t nz,
                               Vec3 const& oxyz, Vec3 const& dxyz)
{
  if (nx == 0 || ny == 0 || nz == 0) return 1;
  // Set origin and spacing
  ox_ = oxyz[0];
  oy_ = oxyz[1];
  oz_ = oxyz[2];
  dx_ = dxyz[0];
  dy_ = dxyz[1];
  dz_ = dxyz[2];
  // Calculate maximum, used when binning
  mx_ = ox_ + ((double)nx * dx_);
  my_ = oy_ + ((double)ny * dy_);
  mz_ = oz_ + ((double)nz * dz_);
  return Allocate3D(nx, ny, nz);
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

// DataSet_3D::Allocate_N_C_D()
int DataSet_3D::Allocate_N_C_D(size_t nx, size_t ny, size_t nz,
                               Vec3 const& cxyz, Vec3 const& dxyz)
{
  // Calculate origin from center coordinates.
  Vec3 oxyz( cxyz[0] + Calc_Origin(nx, dxyz[0]),
             cxyz[1] + Calc_Origin(ny, dxyz[1]),
             cxyz[2] + Calc_Origin(nz, dxyz[2]) );
  return Allocate_N_O_D(nx,ny,nz,oxyz,dxyz);
}

// DataSet_3D::Allocate_X_C_D()
int DataSet_3D::Allocate_X_C_D(Vec3 const& sizes, Vec3 const& center, Vec3 const& dxyz)
{
  // Calculate bin counts
  // TODO: Make size_t
  size_t nx = (size_t)(sizes[0] / dxyz[0]);
  size_t ny = (size_t)(sizes[1] / dxyz[1]);
  size_t nz = (size_t)(sizes[2] / dxyz[2]);
  return Allocate_N_C_D( nx, ny, nz, center, dxyz );
} 
