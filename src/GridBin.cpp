#include "GridBin.h"

/** Set up for orthogonal X-aligned grid with given origin and spacing; calculate maximum. */
void GridBin::Setup_O_D(size_t nx, size_t ny, size_t nz,
                   Vec3 const& oxyzIn, Vec3 const& dxyz)
{
  OXYZ_ = oxyzIn;
  dx_ = dxyz[0];
  dy_ = dxyz[1];
  dz_ = dxyz[2];
  mx_ = OXYZ_[0] + ((double)nx * dx_);
  my_ = OXYZ_[1] + ((double)ny * dy_);
  mz_ = OXYZ_[2] + ((double)nz * dz_);
  nx_ = (double)nx;
  ny_ = (double)ny;
  nz_ = (double)nz;
  voxelvolume_ = dx_ * dy_ * dz_;
  // Set orthogonal unit cell vectors. TODO should these be w.r.t. the offset?
  double ucell[9];
  ucell[0] = (double)nx * dx_;
  ucell[1] = 0;
  ucell[2] = 0;

  ucell[3] = 0;
  ucell[4] = (double)ny * dy_;
  ucell[5] = 0;

  ucell[6] = 0;
  ucell[7] = 0;
  ucell[8] = (double)nz * dz_;

  box_.SetupFromUcell(ucell);
  box_.PrintDebug("GridBin::Setup_O_D");
  SetupInternalPointers();
}

/** Set up for grid with given bins, origin, and box.*/
void GridBin::Setup_O_Box(size_t nxIn, size_t nyIn, size_t nzIn,
                     Vec3 const& oxyzIn, Box const& boxIn)
{
  nx_ = (double)nxIn;
  ny_ = (double)nyIn;
  nz_ = (double)nzIn;
  OXYZ_ = oxyzIn;
  box_ = boxIn;
  box_.PrintDebug("GridBin::Setup_O_Box");
  // Get the 3 individual unit cell vector lengths
  double l_Avec = box_.UnitCell().Row1().Length();
  double l_Bvec = box_.UnitCell().Row2().Length();
  double l_Cvec = box_.UnitCell().Row3().Length();
  // Get spacing from vector length over bins
  dx_ = l_Avec / nx_;
  dy_ = l_Bvec / ny_;
  dz_ = l_Cvec / nz_;
  // Get max from origin plus vector length
  mx_ = OXYZ_[0] + l_Avec;
  my_ = OXYZ_[1] + l_Bvec;
  mz_ = OXYZ_[2] + l_Cvec;
  // Get voxel volume from total grid volume over number of bins.
  voxelvolume_ = boxIn.CellVolume() / (nx_ * ny_ * nz_);
  SetupInternalPointers();
}

