#ifndef INC_DATASET_3D_H
#define INC_DATASET_3D_H
#include "DataSet.h"
#include "CpptrajFile.h"
#include "Vec3.h"
/// Interface for 3D DataSets.
class DataSet_3D : public DataSet {
  public:
    DataSet_3D() {}
    DataSet_3D(DataSet::DataType tIn, int wIn, int pIn) :
      DataSet(tIn, wIn, pIn, 3) {}
    /// Set up grid for given # x, y, and z points.
    virtual int Allocate3D(size_t, size_t, size_t) = 0;
    /// Write 3D data to file.
    virtual void Write3D(CpptrajFile&,int,int,int) const = 0;
    /// \return Data from grid at x/y/z point.
    virtual double GetElement(size_t, size_t, size_t) const = 0;
    /// \return Data from grid.
    virtual double operator[](size_t) const = 0;
    /// \return size of X dimension.
    virtual size_t NX() const = 0;
    /// \return size of Y dimension.
    virtual size_t NY() const = 0;
    /// \return size of Z dimension.
    virtual size_t NZ() const = 0;
    // -------------------------------------------
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
    /// Set up grid origin and spacing (and max for binning).
    int setOriginAndSpacing(size_t,size_t,size_t,double,double,double,double,double,double);
    /// Convert X, Y, and Z coords to index.
    inline bool CalcBins(double,double,double,size_t&,size_t&,size_t&) const;
    inline double DX() const { return dx_; }
    inline double DY() const { return dy_; }
    inline double DZ() const { return dz_; }
    inline double OX() const { return ox_; }
    inline double OY() const { return oy_; }
    inline double OZ() const { return oz_; }
    inline Vec3 BinCorner(size_t,size_t,size_t);
    inline Vec3 BinCenter(size_t,size_t,size_t);
  private:
    double dx_; ///< X grid spacing.
    double dy_; ///< Y grid spacing.
    double dz_; ///< Z grid spacing.
    double ox_; ///< Grid X origin (minimum).
    double oy_; ///< Grid Y origin (minimum).
    double oz_; ///< Grid Z origin (minimum).
    double mx_; ///< Grid X maximum.
    double my_; ///< Grid Y maximum.
    double mz_; ///< Grid Z maximum. 
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
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
// DataSet_3D::CalcBins()
bool DataSet_3D::CalcBins(double x, double y, double z,
                          size_t& i, size_t& j, size_t& k) const
{
  // X
  double xx = x - ox_;
  if (xx >= 0.0 && xx < mx_) {
    // Y
    double yy = y - oy_;
    if (yy >= 0.0 && yy < my_) {
      // Z
      double zz = z - oz_;
      if (zz >= 0.0 && zz < mz_) {
        i = (size_t)(xx / dx_);
        j = (size_t)(yy / dy_);
        k = (size_t)(zz / dz_);
        return true;
      }
    }
  }
  return false;
}
// DataSet_3D::BinCorner()
Vec3 DataSet_3D::BinCorner(size_t i, size_t j, size_t k) {
  return Vec3( (double)i*dx_ + ox_,
               (double)j*dy_ + oy_,
               (double)k*dz_ + oz_ );
}
// DataSet_3D::BinCenter()
Vec3 DataSet_3D::BinCenter(size_t i, size_t j, size_t k) {
  return Vec3( (double)i*dx_ + ox_ + 0.5*dx_,
               (double)j*dy_ + oy_ + 0.5*dy_,
               (double)k*dz_ + oz_ + 0.5*dz_ );
}
#endif
