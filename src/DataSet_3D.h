#ifndef INC_DATASET_3D_H
#define INC_DATASET_3D_H
#include "DataSet.h"
#include "CpptrajFile.h"
#include "Box.h"
#include "GridBin.h"
/// Interface for 3D DataSets.
// FIXME: Use DataSet Dims?
class DataSet_3D : public DataSet {
  public:
    DataSet_3D() : gridBin_(0) {}
    virtual ~DataSet_3D(); // Virtual since this class is inherited.
    DataSet_3D(DataSet_3D const&);
    DataSet_3D& operator=(DataSet_3D const&);
    DataSet_3D(DataSet::DataType tIn, TextFormat const& fIn) :
      DataSet(tIn, GRID_3D, fIn, 3), gridBin_(0) {}
    // TODO enable append?
    int Append(DataSet*) { return 1; }
    int Allocate(SizeArray const&) { return 1; } // TODO enable?
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
    /// \return grid index
    virtual long int CalcIndex(size_t, size_t, size_t) const = 0;
    /// Calculate bins for given grid index
    virtual void ReverseIndex(long int, size_t&, size_t&, size_t&) const = 0;
    /// Increment specified voxel by given amount.
    virtual void UpdateVoxel(long int, double) = 0;
    // -------------------------------------------
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
    /*
    double Coord(unsigned int d, size_t p) const {
      long int idx[3];
      for (unsigned int i = 0; i < 3; i++)
        if (i == d) idx[i] = p; else idx[i] = 0;
      Vec3 crd = gridBin_->Corner(idx[0], idx[1], idx[2]);
      return crd[d];
    }*/
    // -------------------------------------------
    /// Set up grid from dims, origin, and spacing.
    int Allocate_N_O_D(size_t,size_t,size_t,Vec3 const&,Vec3 const&);
    /// Set up grid from dims, center, and spacing.
    int Allocate_N_C_D(size_t,size_t,size_t,Vec3 const&,Vec3 const&);
    /// Set up grid from sizes, center, and spacing.
    int Allocate_X_C_D(Vec3 const&,Vec3 const&,Vec3 const&);
    /// Set up grid from dims, origin, and box.
    int Allocate_N_O_Box(size_t,size_t,size_t, Vec3 const&, Box const&);
    /// Print grid info.
    void GridInfo() const;
    // -------------------------------------------
    GridBin const& Bin() const { return *gridBin_; }
  private:
    /// Check if grid dimension is even; if not, increment it by 1.
    static void CheckEven(size_t&, char);
    /// Set up grid for given # x, y, and z points.
    // TODO: Make public if grids will be used for other than binning.
    virtual int Allocate3D(size_t, size_t, size_t) = 0;

    GridBin* gridBin_; ///< Used to calculate bins/coords depending on grid type.
};
#endif
