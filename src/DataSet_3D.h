#ifndef INC_DATASET_3D_H
#define INC_DATASET_3D_H
#include "DataSet.h"
#include "CpptrajFile.h"
#include "Vec3.h"
/// Interface for 3D DataSets.
// FIXME: Use DataSet Dims?
class DataSet_3D : public DataSet {
  public:
    DataSet_3D() {}
    DataSet_3D(DataSet::DataType tIn, int wIn, int pIn) :
      DataSet(tIn, wIn, pIn, 3) {}
    /// Write 3D data to file.
    virtual void Write3D(CpptrajFile&,int,int,int) const = 0;
    /// \return Data from grid at x/y/z point.
    virtual double GetElement(int, int, int) const = 0;
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
    /// Set up grid from dims, origin, and spacing.
    int Allocate_N_O_D(size_t,size_t,size_t,Vec3 const&,Vec3 const&);
    /// Set up grid from dims, center, and spacing.
    int Allocate_N_C_D(size_t,size_t,size_t,Vec3 const&,Vec3 const&);
    /// Set up grid from sizes, center, and spacing.
    int Allocate_X_C_D(Vec3 const&,Vec3 const&,Vec3 const&);
    /// Convert X, Y, and Z coords to index.
    inline bool CalcBins(double,double,double,int&,int&,int&) const;
    inline double DX() const { return dx_; }
    inline double DY() const { return dy_; }
    inline double DZ() const { return dz_; }
    inline double OX() const { return ox_; }
    inline double OY() const { return oy_; }
    inline double OZ() const { return oz_; }
    inline double MX() const { return mx_; }
    inline double MY() const { return my_; }
    inline double MZ() const { return mz_; }
    inline Vec3 BinCorner(int,int,int);
    inline Vec3 BinCenter(int,int,int);
  private:
    /// Check if grid dimension is even; if not, increment it by 1.
    static void CheckEven(size_t&, char);
    /// Set up grid for given # x, y, and z points.
    // TODO: Make public if grids will be used for other than binning.
    virtual int Allocate3D(size_t, size_t, size_t) = 0;

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
// DataSet_3D::CalcBins()
bool DataSet_3D::CalcBins(double x, double y, double z,
                          int& i, int& j, int& k) const
{
  // X
  if (x >= ox_ && x < mx_) {
    // Y
    if (y >= oy_ && y < my_) {
      // Z
      if (z >= oz_ && z < mz_) {
        i = (int)((x-ox_) / dx_);
        j = (int)((y-oy_) / dy_);
        k = (int)((z-oz_) / dz_);
        return true;
      }
    }
  }
  return false;
}
// DataSet_3D::BinCorner()
Vec3 DataSet_3D::BinCorner(int i, int j, int k) {
  return Vec3( (double)i*dx_ + ox_,
               (double)j*dy_ + oy_,
               (double)k*dz_ + oz_ );
}
// DataSet_3D::BinCenter()
Vec3 DataSet_3D::BinCenter(int i, int j, int k) {
  return Vec3( (double)i*dx_ + ox_ + 0.5*dx_,
               (double)j*dy_ + oy_ + 0.5*dy_,
               (double)k*dz_ + oz_ + 0.5*dz_ );
}
#endif
