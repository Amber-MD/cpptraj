#ifndef INC_GRIDBIN_H
#define INC_GRIDBIN_H
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
    bool Calc(double, double, double, size_t&, size_t&, size_t&) const;

    /// Given coordinates, set corresponding bin indices; no bounds check.
    virtual void Indices(double, double, double, long int&, long int&, long int&) const = 0;
    /// \return coordinates of bin for given indices; no bound check.
    virtual Vec3 Corner(long int, long int, long int) const = 0;
    /// \return coordinates of bin center for given indices; no bounds check.
    virtual Vec3 Center(long int, long int, long int) const = 0;
    /// \return unit cell matrix. // TODO: Make const&?
    virtual Matrix_3x3 Ucell() const = 0;
    /// \return true if GridBin type is orthogonal.
    virtual bool IsOrthoGrid() const = 0;
    /// \return Voxel volume.
    virtual double VoxelVolume() const = 0;
    /// \return a copy of this GridBin.
    virtual GridBin* Copy() const = 0;
    /// \return Grid origin.
    Vec3 const& GridOrigin() const { return OXYZ_; }
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

#endif
