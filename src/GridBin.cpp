#include "GridBin.h"
#include <cmath> // ceil

/** Set voxel volume from total grid volume over number of bins. */
void GridBin::set_voxel_volume() {
  voxelvolume_ = box_.CellVolume() / (nx_ * ny_ * nz_);
}

/** Put given grid sizes into an array. */
static inline GridBin::SizeArray getGridSizes(size_t nx, size_t ny, size_t nz)
{
  GridBin::SizeArray gridSizes(3);
  gridSizes[0] = nx;
  gridSizes[1] = ny;
  gridSizes[2] = nz;
  return gridSizes;
}

/** Set up for grid with given bins, origin, and box.*/
GridBin::SizeArray GridBin::Setup_Sizes_Origin_Box(size_t nxIn, size_t nyIn, size_t nzIn,
                                                   Vec3 const& oxyzIn, Box const& boxIn)
{
  // Set grid dimensions
  nx_ = (double)nxIn;
  ny_ = (double)nyIn;
  nz_ = (double)nzIn;
  OXYZ_ = oxyzIn;
  // Set grid box and internal pointers based on box.
  box_ = boxIn;
  box_.PrintDebug("GridBin::Setup_Sizes_Origin_Box");
  SetupInternalPointers();
  set_voxel_volume();
  // Get the 3 individual unit cell vector lengths
  double l_Avec = box_.UnitCell().Row1().Length();
  double l_Bvec = box_.UnitCell().Row2().Length();
  double l_Cvec = box_.UnitCell().Row3().Length();
  // Get spacing from vector length over bins
  dx_ = l_Avec / nx_;
  dy_ = l_Bvec / ny_;
  dz_ = l_Cvec / nz_;
  // Set origin and max
  SetOrigin( oxyzIn );

  return getGridSizes(nxIn, nyIn, nzIn);
}

/** Set up for orthogonal X-aligned grid with given origin and spacing; calculate maximum. */
GridBin::SizeArray GridBin::Setup_Sizes_Origin_Spacing(size_t nx, size_t ny, size_t nz,
                                                       Vec3 const& oxyzIn, Vec3 const& dxyz)
{
  // Set grid dimensions
  nx_ = (double)nx;
  ny_ = (double)ny;
  nz_ = (double)nz;
  // Set grid spacings
  dx_ = dxyz[0];
  dy_ = dxyz[1];
  dz_ = dxyz[2];
  // Set grid box and internal pointers based on box.
  box_.SetupFromXyzAbg( nx_ * dx_, ny_ * dy_, nz_ * dz_, 90.0, 90.0, 90.0 );
  box_.PrintDebug("GridBin::Setup_Sizes_Origin_Spacing");
  SetupInternalPointers();
  set_voxel_volume();
  // Set origin and max
  SetOrigin( oxyzIn );

  return getGridSizes(nx, ny, nz);
}

/** Set up for orthogonal X-aligned grid with given center and spacing; calculate maximum. */
GridBin::SizeArray GridBin::Setup_Sizes_Center_Spacing(size_t nx, size_t ny, size_t nz,
                                                       Vec3 const& cxyzIn, Vec3 const& dxyz)
{
  // Set grid dimensions
  nx_ = (double)nx;
  ny_ = (double)ny;
  nz_ = (double)nz;
  // Set grid spacings
  dx_ = dxyz[0];
  dy_ = dxyz[1];
  dz_ = dxyz[2];
  // Set grid box and internal pointers based on box.
  box_.SetupFromXyzAbg( nx_ * dx_, ny_ * dy_, nz_ * dz_, 90.0, 90.0, 90.0 );
  box_.PrintDebug("GridBin::Setup_Sizes_Center_Spacing");
  SetupInternalPointers();
  set_voxel_volume();
  // Set origin and max
  SetOriginFromCenter( cxyzIn );

  return getGridSizes(nx, ny, nz);
}

/** Set up for orthogonal X-aligned grid with given lengths, center and spacing. */
GridBin::SizeArray GridBin::Setup_Lengths_Center_Spacing(Vec3 const& lengths, Vec3 const& center,
                                                         Vec3 const& dxyz)
{
  // Calculate bin counts
  size_t nx = (size_t)ceil(lengths[0] / dxyz[0]);
  size_t ny = (size_t)ceil(lengths[1] / dxyz[1]);
  size_t nz = (size_t)ceil(lengths[2] / dxyz[2]);
  return Setup_Sizes_Center_Spacing(nx, ny, nz, center, dxyz);
}
