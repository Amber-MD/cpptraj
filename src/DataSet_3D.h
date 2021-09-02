#ifndef INC_DATASET_3D_H
#define INC_DATASET_3D_H
#include <cstddef> // size_t
#include "DataSet.h"
#include "GridBin.h"
class Box;
/// Interface for 3D DataSets.
/** Compile with -DDEBUG_DATASET_3D to print debug info during allocation. */
// TODO: Use DataSet Dims?
class DataSet_3D : public DataSet {
  public:
    /// CONSTRUCTOR
    DataSet_3D();
    virtual ~DataSet_3D(); // Virtual since this class is inherited.
    /// COPY CONSTRUCTOR
    DataSet_3D(DataSet_3D const&);
    /// ASSIGNMENT
    DataSet_3D& operator=(DataSet_3D const&);
    /// CONSTRUCTOR - type, format
    DataSet_3D(DataSet::DataType tIn, TextFormat const& fIn) :
      DataSet(tIn, GRID_3D, fIn, 3) {}
    // ----- DataSet -----------------------------
    // TODO enable append?
    int Append(DataSet*) { return 1; }
    int Allocate(SizeArray const&) { return 1; } // TODO enable?
    // TODO: Remove this. Only needed by DataSet.h
    void Add(size_t,const void*) { }
#   ifdef MPI
    /// Sum grid across ranks to master, ensure orientation from final rank is sent to master.
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    // -------------------------------------------
    /// \return Data from grid at x/y/z point.
    virtual double GetElement(size_t, size_t, size_t) const = 0;
    /// Set grid to value
    virtual void SetGrid(size_t, double) = 0;
    /// Increment the specified grid point by value
    virtual void IncrementElement(size_t, size_t, size_t, double) = 0;
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
    /// Divide all elements by the given scalar
    virtual void operator/=(double) = 0;
    // -------------------------------------------
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

    /// Move grid center
    void SetGridCenter(Vec3 const& cxyz) { gridBin_.SetOriginFromCenter( cxyz); }
    /// Set the grid unit cell
    void Assign_Grid_UnitCell(Matrix_3x3 const& ucell) { gridBin_.Assign_UnitCell(ucell); }
    /// Rotate the grid using the given rotation matrix.
    void Rotate_3D_Grid(Matrix_3x3 const& Rot) { gridBin_.RotateGrid(Rot); }
    /// X-align the grid
    void Xalign_3D_Grid() { gridBin_.X_align_grid(); }
    /// Print grid info.
    void GridInfo() const;
    /// \return GridBin interface
    GridBin const& Bin() const { return gridBin_; }
    // -------------------------------------------
  protected:
#   ifdef MPI
    /** Used by inheriting class to sum grid across ranks to master. */
    virtual int SyncGrid(size_t, std::vector<int> const&, Parallel::Comm const&) = 0;
#   endif
  private:
    /// Check if grid dimension is even; if not, increment it by 1.
    //static void CheckEven(size_t&, char);
    /// Set up grid for given # x, y, and z points.
    // TODO: Make public if grids will be used for other than binning.
    virtual int Allocate3D(size_t, size_t, size_t) = 0;
    /// \return Origin coords from center, spacing, and sizes
    //static Vec3 calcOriginFromCenter(Vec3 const&, double, double, double, size_t, size_t, size_t);

    GridBin gridBin_; ///< Used to calculate bins/coords depending on grid orientation.
};
#endif
