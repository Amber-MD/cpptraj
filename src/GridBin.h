#ifndef INC_GRIDBIN_H
#define INC_GRIDBIN_H
#include <cstddef> // size_t
#include <vector>
#include "Box.h"
/// Class used to describe 3D grid spatial orientation, as well as perform binning on/get voxel coords.
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
    /// \return true if Grid is X-aligned
    inline bool IsXalignedGrid() const { return box_.Is_X_Aligned(); }
    /// \return Voxel volume.
    inline double VoxelVolume() const { return voxelvolume_; }
    /// \return a copy of this GridBin.
    //virtual GridBin* Copy() const = 0;
    /// \return Grid origin coordinates.
    inline Vec3 const& GridOrigin() const { return OXYZ_; }
    /// \return Grid center coordinates.
    inline Vec3 GridCenter() const;

    // TODO are these spacing routines needed?
    inline double DX() const { return dx_; }
    inline double DY() const { return dy_; }
    inline double DZ() const { return dz_; }

    /// Print debug info
    void PrintDebug(const char*) const;

    // Set up routines.
    /// Set grid origin
    inline void SetOrigin(Vec3 const&);
    /// Set grid origin via given center point.
    inline void SetOriginFromCenter(Vec3 const&);
    /// Apply rotation matrix to grid unit cell vectors.
    void RotateGrid(Matrix_3x3 const&);
    /// Make the grid X-aligned
    void X_align_grid();
    /// Type for returning grid dimensions in terms of bins from setup routines
    typedef std::vector<size_t> SizeArray;
    /// Set up for grid with given bins, origin, and box.
    SizeArray Setup_Sizes_Origin_Box(size_t, size_t, size_t, Vec3 const&, Box const&);
    /// Set up for orthogonal X-aligned grid with given bins, origin and spacing; calculate maximum.
    SizeArray Setup_Sizes_Origin_Spacing(size_t, size_t, size_t, Vec3 const&, Vec3 const&);
    /// Set up for orthogonal X-aligned grid with given bins, center and spacing; calculate maximum.
    SizeArray Setup_Sizes_Center_Spacing(size_t, size_t, size_t, Vec3 const&, Vec3 const&);
    /// Set up for orthogonal X-aligned grid with given lengths, center, and spacing.
    SizeArray Setup_Lengths_Center_Spacing(Vec3 const&, Vec3 const&, Vec3 const&);
    /// Assign new unit cell vectors to grid
    void Assign_UnitCell( Matrix_3x3 const& );
#   ifdef MPI
    /// Assume split across trajectory; ensure master has orientation from final rank
    int Sync(Parallel::Comm const&);
#   endif
  private:
    inline bool Calc_ortho(double, double, double, size_t&, size_t&, size_t&) const;
    inline bool Calc_nonortho(double, double, double, size_t&, size_t&, size_t&) const;
    inline void Indices_ortho(double, double, double, long int&, long int&, long int&) const;
    inline void Indices_nonortho(double, double, double, long int&, long int&, long int&) const;
    inline Vec3 Corner_ortho(long int, long int, long int) const;
    inline Vec3 Corner_nonortho(long int, long int, long int) const;
    inline Vec3 Center_ortho(long int, long int, long int) const;
    inline Vec3 Center_nonortho(long int, long int, long int) const;
    inline void SetupInternalPointers();
    void set_voxel_volume();

    Vec3 OXYZ_;           ///< Grid origin in Cartesian coordinates.
    double dx_, dy_, dz_; ///< Grid spacing for binning (Ang., Cartesian).
    double mx_, my_, mz_; ///< Grid max for binning (Ang., Cartesian, orthogonal).
    double nx_, ny_, nz_; ///< Number of bins in double precision (for nonortho binning).
    double voxelvolume_;  ///< Volume of a single voxel (Ang^3).
    Box box_;             ///< Contain grid unit cell vectors, frac. vectors, volume.
    /// Internal pointer to the correct Calc routine
    bool (GridBin::*CalcPtr_)(double, double, double, size_t&, size_t&, size_t&) const;
    /// Internal pointer to the correct Indices routine
    void (GridBin::*IndicesPtr_)(double, double, double, long int&, long int&, long int&) const;
    /// Internal pointer to the correct Corner routine
    Vec3 (GridBin::*CornerPtr_)(long int, long int, long int) const;
    /// Internal pointer to the correct Center routine
    Vec3 (GridBin::*CenterPtr_)(long int, long int, long int) const;
};
// -----------------------------------------------------------------------------

/** Interface to the Calc routine appropriate for the grid type. */
bool GridBin::Calc(double x, double y, double z, size_t& i, size_t& j, size_t& k) const {
  return ((*this).*(CalcPtr_))(x, y, z, i, j, k);
}

/** Interface to the Indices routine appropriate for the grid type. */
void GridBin::Indices(double x, double y, double z, long int& i, long int& j, long int& k) const
{
  ((*this).*(IndicesPtr_))(x, y, z, i, j, k);
}

/** Interface to the Corner routine appropriate for the grid type. */
Vec3 GridBin::Corner(long int i, long int j, long int k) const
{
  return ((*this).*(CornerPtr_))(i, j, k);
}

/** Interface to the Center routine appropriate for the grid type. */
Vec3 GridBin::Center(long int i, long int j, long int k) const
{
  return ((*this).*(CenterPtr_))(i, j, k);
}

/** \return true if given coordinates are on grid; set corresponding bin indices. */
bool GridBin::Calc_ortho(double x, double y, double z, size_t& i, size_t& j, size_t& k) const
{
  //if (box_.Is_X_Aligned_Ortho()) {
    //mprintf("DEBUG: X-aligned Calc\n");
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
  //} else {
  return false;
}

/** \return true if given coordinates are on grid; set corresponding bin indices. */
bool GridBin::Calc_nonortho(double x, double y, double z, size_t& i, size_t& j, size_t& k) const
{
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
  //}
  return false;
}

/** Given coordinates, set corresponding bin indices; no bounds check. */
void GridBin::Indices_ortho(double x, double y, double z, long int& i, long int& j, long int& k) const
{
  //if (box_.Is_X_Aligned_Ortho()) {
    // X-aligned and orthogonal
    i = (long int)((x-OXYZ_[0]) / dx_);
    j = (long int)((y-OXYZ_[1]) / dy_);
    k = (long int)((z-OXYZ_[2]) / dz_);
  //} else {
}

/** Given coordinates, set corresponding bin indices; no bounds check. */
void GridBin::Indices_nonortho(double x, double y, double z, long int& i, long int& j, long int& k) const
{
    // Not X-aligned or non-orthogonal
    Vec3 frac = box_.FracCell() * Vec3(x - OXYZ_[0], y - OXYZ_[1], z - OXYZ_[2]);
    i = (long int)(frac[0] * nx_);
    j = (long int)(frac[1] * ny_);
    k = (long int)(frac[2] * nz_);
  //}
}


/** \return coordinates of bin for given indices; no bound check. */
Vec3 GridBin::Corner_ortho(long int i, long int j, long int k) const
{
  //if (box_.Is_X_Aligned_Ortho()) {
    // X-aligned and orthogonal
    return Vec3((double)i*dx_+OXYZ_[0],
                (double)j*dy_+OXYZ_[1],
                (double)k*dz_+OXYZ_[2]);
  //} else {
}

/** \return coordinates of bin for given indices; no bound check. */
Vec3 GridBin::Corner_nonortho(long int i, long int j, long int k) const
{
    // Not X-aligned or non-orthogonal
    Vec3 frac( (double)i / nx_, (double)j / ny_, (double)k / nz_ );
    return box_.UnitCell().TransposeMult( frac ) + OXYZ_;
  //}
}

/** \return coordinates of bin center for given indices; no bounds check. */
Vec3 GridBin::Center_ortho(long int i, long int j, long int k) const
{
  //if (box_.Is_X_Aligned_Ortho()) {
    // X-aligned and orthogonal
    return Vec3((double)i*dx_+OXYZ_[0]+0.5*dx_,
                (double)j*dy_+OXYZ_[1]+0.5*dy_,
                (double)k*dz_+OXYZ_[2]+0.5*dz_);
  //} else {
}

/** \return coordinates of bin center for given indices; no bounds check. */
Vec3 GridBin::Center_nonortho(long int i, long int j, long int k) const
{
    // Not X-aligned or non-orthogonal
    Vec3 frac_half((1.0 + 2.0 * (double)i) / (2.0 * nx_),  //(0.5 * (1.0 / nx_)) + ((double)i / nx_),
                   (1.0 + 2.0 * (double)j) / (2.0 * ny_), 
                   (1.0 + 2.0 * (double)k) / (2.0 * nz_));
    return box_.UnitCell().TransposeMult( frac_half ) + OXYZ_;
  //}
}

/** Set up function internal pointers. */
void GridBin::SetupInternalPointers() {
  if (box_.Is_X_Aligned_Ortho()) {
    CalcPtr_ = &GridBin::Calc_ortho;
    IndicesPtr_ = &GridBin::Indices_ortho;
    CornerPtr_ = &GridBin::Corner_ortho;
    CenterPtr_ = &GridBin::Center_ortho;
  } else {
    CalcPtr_ = &GridBin::Calc_nonortho;
    IndicesPtr_ = &GridBin::Indices_nonortho;
    CornerPtr_ = &GridBin::Corner_nonortho;
    CenterPtr_ = &GridBin::Center_nonortho;
  }
}

/** Set new origin for grid, update max. */
void GridBin::SetOrigin(Vec3 const& newOxyz) {
  OXYZ_ = newOxyz;
  if (box_.Is_X_Aligned_Ortho()) {
    mx_ = OXYZ_[0] + (nx_ * dx_);
    my_ = OXYZ_[1] + (ny_ * dy_);
    mz_ = OXYZ_[2] + (nz_ * dz_);
  } else {
    // Get vector pointing from coord origin to max grid
    // TODO is this necessary? mx etc not used unless X aligned...
    Vec3 maxVec = box_.UnitCell().TransposeMult( Vec3(1.0, 1.0, 1.0) );
    mx_ = OXYZ_[0] + maxVec[0];
    my_ = OXYZ_[1] + maxVec[1];
    mz_ = OXYZ_[2] + maxVec[2];
  }
}

/** Set new origin for grid based on location of center, update max. */
void GridBin::SetOriginFromCenter(Vec3 const& newCxyz) {
  // Get vector pointing from coord origin to center
  Vec3 centerVec;
  if (box_.Is_X_Aligned_Ortho()) {
    centerVec = Vec3(box_.Param(Box::X)/2.0, box_.Param(Box::Y)/2.0, box_.Param(Box::Z)/2.0);
  } else {
    centerVec = box_.UnitCell().TransposeMult( Vec3(0.5, 0.5, 0.5) );
  }
  OXYZ_ = newCxyz - centerVec;
  mx_ = newCxyz[0] + centerVec[0];
  my_ = newCxyz[1] + centerVec[1];
  mz_ = newCxyz[2] + centerVec[2];
}

/** \return Grid center coordinates. */
Vec3 GridBin::GridCenter() const {
  // Get vector pointing from coord origin to center
  Vec3 centerVec;
  if (box_.Is_X_Aligned_Ortho()) {
    centerVec = Vec3(box_.Param(Box::X)/2.0, box_.Param(Box::Y)/2.0, box_.Param(Box::Z)/2.0);
  } else {
    centerVec = box_.UnitCell().TransposeMult( Vec3(0.5, 0.5, 0.5) );
  }
  return OXYZ_ + centerVec;
}
#endif
