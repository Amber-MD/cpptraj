#ifndef INC_GRIDBIN_H
#define INC_GRIDBIN_H
#include <cstddef> // size_t
#include "Box.h"
/// Class used to perform binning on/get voxel coords of 3D grids.
class GridBin {
  public:
    /// CONSTRUCTOR
    GridBin() : OXYZ_(0.0),
                dx_(-1.0), dy_(-1.0), dz_(-1.0), mx_(0),  my_(0), mz_(0),
                nx_(0), ny_(0), nz_(0), voxelvolume_(0)
                {}
    /// \return true if given coordinates are on grid; set corresponding bin indices.
    inline bool Calc(double, double, double, size_t&, size_t&, size_t&) const;
    /// Given coordinates, set corresponding bin indices; no bounds check.
    inline void Indices(double, double, double, long int&, long int&, long int&) const;
    /// \return coordinates of bin for given indices; no bound check.
    inline Vec3 Corner(long int, long int, long int) const;
    /// \return coordinates of bin center for given indices; no bounds check.
    inline Vec3 Center(long int, long int, long int) const;
    /// \return Grid unit cell matrix. TODO just use GridBox?
    inline Matrix_3x3 const& Ucell() const { return box_.UnitCell(); }
    /// \return Grid box
    inline Box const& GridBox() const { return box_; }
    /// \return true if Grid is X-aligned and orthogonal. TODO this may not be necessary
    inline bool IsOrthoGrid() const { return box_.Is_X_Aligned_Ortho(); }
    /// \return Voxel volume.
    inline double VoxelVolume() const { return voxelvolume_; }
    /// \return a copy of this GridBin.
    //virtual GridBin* Copy() const = 0;
    /// \return Grid origin.
    inline Vec3 const& GridOrigin() const { return OXYZ_; }

    // TODO are these spacing routines needed?
    inline double DX() const { return dx_; }
    inline double DY() const { return dy_; }
    inline double DZ() const { return dz_; }

    // Set up routines.
    /// Set up for orthogonal X-aligned grid with given origin and spacing; calculate maximum.
    inline void Setup_O_D(size_t, size_t, size_t, Vec3 const&, Vec3 const&);
    /// Set up for grid with given bins, origin, and box.
    inline void Setup_O_Box(size_t, size_t, size_t, Vec3 const&, Box const&);
  protected:
    Vec3 OXYZ_;           ///< Grid origin.
    double dx_, dy_, dz_; ///< Grid spacing (Ang., Cartesian, orthogonal).
    double mx_, my_, mz_; ///< Grid max (Ang., Cartesian, orthogonal).
    double nx_, ny_, nz_; ///< Number of bins in double precision (nonortho).
    double voxelvolume_;  ///< Volume of a single voxel (Ang^3).
    Box box_;             ///< Contain grid unit cell vectors, frac. vectors, volume.
};
// -----------------------------------------------------------------------------

/** \return true if given coordinates are on grid; set corresponding bin indices. */
bool GridBin::Calc(double x, double y, double z, size_t& i, size_t& j, size_t& k) const
{
  if (box_.Is_X_Aligned_Ortho()) {
    // X-aligned and orthogonal
    if (x >= OXYZ_[0] && x < mx_) { // X
      if (y >= OXYZ_[1] && y < my_) { // Y
        if (z >= OXYZ_[2] && z < mz_) { // Z
          i = (size_t)((x-OXYZ_[0]) / dx_);
          j = (size_t)((y-OXYZ_[1]) / dy_);
          k = (size_t)((z-OXYZ_[2]) / dz_);
          return true;
        }
      }
    }
  } else {
    // Not X-aligned or non-orthogonal
    Vec3 frac = box_.FracCell() * Vec3(x - OXYZ_[0], y - OXYZ_[1], z - OXYZ_[2]);
    if (frac[0] >= 0.0 && frac[0] < 1.0) {
      if (frac[1] >= 0.0 && frac[1] < 1.0) {
        if (frac[2] >= 0.0 && frac[2] < 1.0) {
          i = (size_t)(frac[0] * nx_);
          j = (size_t)(frac[1] * ny_);
          k = (size_t)(frac[2] * nz_);
          return true;
        }
      }
    }
  }
  return false;
}

/** Given coordinates, set corresponding bin indices; no bounds check. */
void GridBin::Indices(double x, double y, double z, long int& i, long int& j, long int& k) const
{
  if (box_.Is_X_Aligned_Ortho()) {
    // X-aligned and orthogonal
    i = (long int)((x-OXYZ_[0]) / dx_);
    j = (long int)((y-OXYZ_[1]) / dy_);
    k = (long int)((z-OXYZ_[2]) / dz_);
  } else {
    // Not X-aligned or non-orthogonal
    Vec3 frac = box_.FracCell() * Vec3(x - OXYZ_[0], y - OXYZ_[1], z - OXYZ_[2]);
    i = (long int)(frac[0] * nx_);
    j = (long int)(frac[1] * ny_);
    k = (long int)(frac[2] * nz_);
  }
}


/** \return coordinates of bin for given indices; no bound check. */
Vec3 GridBin::Corner(long int i, long int j, long int k) const
{
  if (box_.Is_X_Aligned_Ortho()) {
    // X-aligned and orthogonal
    return Vec3((double)i*dx_+OXYZ_[0],
                (double)j*dy_+OXYZ_[1],
                (double)k*dz_+OXYZ_[2]);
  } else {
    // Not X-aligned or non-orthogonal
    Vec3 frac( (double)i / nx_, (double)j / ny_, (double)k / nz_ );
    return box_.UnitCell().TransposeMult( frac ) + OXYZ_;
  }
}

/** \return coordinates of bin center for given indices; no bounds check. */
Vec3 GridBin::Center(long int i, long int j, long int k) const
{
  if (box_.Is_X_Aligned_Ortho()) {
    // X-aligned and orthogonal
    return Vec3((double)i*dx_+OXYZ_[0]+0.5*dx_,
                (double)j*dy_+OXYZ_[1]+0.5*dy_,
                (double)k*dz_+OXYZ_[2]+0.5*dz_);
  } else {
    // Not X-aligned or non-orthogonal
    Vec3 frac_half((1.0 + 2.0 * (double)i) / (2.0 * nx_),  //(0.5 * (1.0 / nx_)) + ((double)i / nx_),
                   (1.0 + 2.0 * (double)j) / (2.0 * ny_), 
                   (1.0 + 2.0 * (double)k) / (2.0 * nz_));
    return box_.UnitCell().TransposeMult( frac_half ) + OXYZ_;
  }
}

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
}

#endif
