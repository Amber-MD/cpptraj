#include "GridBin.h"
#include "CpptrajStdio.h"
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

/** Apply given unit cell vectors to this grid, keeping the current center
  * location. Internal pointers are reset for the new orientation, and the
  * potentially new volume and origin coordinates are calculated.
  */
void GridBin::Assign_UnitCell( Matrix_3x3 const& unitCell ) {
  // Save the grid center coords TODO should this just always be saved? Can we skip?
  Vec3 gridCtrXyz = GridCenter();

  box_.SetupFromUcell( unitCell.Dptr() );
  SetupInternalPointers();
  set_voxel_volume();

  // Set origin and max
  SetOriginFromCenter( gridCtrXyz );
}

/** Print debug info. */
void GridBin::PrintDebug(const char* title) const
{
  box_.PrintDebug(title);
  mprintf("DEBUG: %s: origin xyz : %12.4f %12.4f %12.4f\n", title, OXYZ_[0], OXYZ_[1], OXYZ_[2]);
  mprintf("DEBUG: %s: spacings   : %12.4f %12.4f %12.4f\n", title, dx_, dy_, dz_);
  mprintf("DEBUG: %s: max xyz    : %12.4f %12.4f %12.4f\n", title, mz_, my_, mz_);
  mprintf("DEBUG: %s: sizes      : %12.0f %12.0f %12.0f\n", title, nx_, ny_, nz_);
  mprintf("DEBUG: %s: voxel vol. : %12.4f\n", title, voxelvolume_);
}

/// Apply translation, rotation, translation
/*
static inline void T_R_T(double* point, Vec3 const& t1, Matrix_3x3 const& R, Vec3 const& t2)
{
  double x = point[0] + t1[0];
  double y = point[1] + t1[1];
  double z = point[2] + t1[2];
  point[0] = x*R[0] + y*R[1] + z*R[2] + t2[0];
  point[1] = x*R[3] + y*R[4] + z*R[5] + t2[1];
  point[2] = x*R[6] + y*R[7] + z*R[8] + t2[2];
}
*/

/** Apply rotation to grid unit cell vectors. */
void GridBin::RotateGrid(Matrix_3x3 const& Rot)
{
  // Save the grid center coords TODO should this just always be saved?
  Vec3 gridCtrXyz = GridCenter();
  //mprintf("DEBUG: Original oxyz= %f %f %f\n", OXYZ_[0], OXYZ_[1], OXYZ_[2]);
  // Rotate grid unit cell
  box_.RotateUcell(Rot);
  // Update internal pointers based on new cell orientation
  SetupInternalPointers();
  //set_voxel_volume(); // Voxel volume should be unchanged
  // Update origin by setting the grid back at the original center 
  // with the new grid unit cell.
  SetOriginFromCenter( gridCtrXyz );
  //mprintf("DEBUG: New oxyz= %f %f %f\n", OXYZ_[0], OXYZ_[1], OXYZ_[2]);
}

/** X-align the current grid. */
void GridBin::X_align_grid() {
  // Save the grid center coords
  Vec3 gridCtrXyz = GridCenter();
  // Create X-aligned box based on current xyz abg.
  double newbox[6];
  newbox[Box::X] = box_.Param(Box::X);
  newbox[Box::Y] = box_.Param(Box::Y);
  newbox[Box::Z] = box_.Param(Box::Z);
  newbox[Box::ALPHA] = box_.Param(Box::ALPHA);
  newbox[Box::BETA] = box_.Param(Box::BETA);
  newbox[Box::GAMMA] = box_.Param(Box::GAMMA);
  box_.AssignFromXyzAbg( newbox );
  // Update internal pointers based on new cell orientation
  SetupInternalPointers();
  // Update origin by setting the grid back at the original center 
  // with the new grid unit cell.
  SetOriginFromCenter( gridCtrXyz );
}

#ifdef MPI
/** Assuming the dataset was split across given comm and that the
  * final rank has the final orientation, ensure master rank
  * has that orientation/location. Should already have same
  * spacing, bin sizes, and volume.
  */
int GridBin::Sync(Parallel::Comm const& commIn) {
  // If only 1 rank no need for this.
  if (commIn.Size() < 2) return 0;
  // Determine final rank
  int finalRank = commIn.Size() - 1;
  if (commIn.Master()) {
    commIn.Recv( OXYZ_.Dptr(), 3, MPI_DOUBLE, finalRank, 2000 );
    commIn.Recv( &mx_,         1, MPI_DOUBLE, finalRank, 2001 );
    commIn.Recv( &my_,         1, MPI_DOUBLE, finalRank, 2002 );
    commIn.Recv( &mz_,         1, MPI_DOUBLE, finalRank, 2003 );
    box_.RecvBox( finalRank, commIn );
    // Update internal pointers based on new cell orientation
    SetupInternalPointers();
  } else if (commIn.Rank() == finalRank) {
    commIn.Send( OXYZ_.Dptr(), 3, MPI_DOUBLE, 0,         2000 );
    commIn.Send( &mx_,         1, MPI_DOUBLE, 0,         2001 );
    commIn.Send( &my_,         1, MPI_DOUBLE, 0,         2002 );
    commIn.Send( &mz_,         1, MPI_DOUBLE, 0,         2003 );
    box_.SendBox( 0, commIn );
  }
  commIn.Barrier();
  return 0;
}
#endif
