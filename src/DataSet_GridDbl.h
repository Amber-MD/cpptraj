#ifndef INC_DATASET_GRDDBL_H
#define INC_DATASET_GRDDBL_H
#include "DataSet_3D.h"
#include "Grid.h"
/// Double-precision three-dimensional grid.
class DataSet_GridDbl : public DataSet_3D {
  public:
    DataSet_GridDbl() : DataSet_3D(GRID_FLT, TextFormat(TextFormat::DOUBLE, 12, 4)) {}
    DataSet_GridDbl(DataSet_GridDbl const&);
    double& operator[](size_t idx)              { return grid_[idx];          }
    static DataSet* Alloc()       { return (DataSet*)new DataSet_GridDbl();  }
    Grid<double> const& InternalGrid()    const { return grid_; }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return grid_.size();        }
#   ifdef MPI
    // FIXME: Currently just sums up. Should this be a separate Sync function?
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info()                          const { return;                     }
    void WriteBuffer(CpptrajFile&,SizeArray const&) const;
    // ----- DataSet_3D functions ----------------
    int Allocate3D(size_t x,size_t y,size_t z) { return grid_.resize(x,y,z); }
    double GetElement(int x,int y,int z) const { return grid_.element(x,y,z); }
    double operator[](size_t idx)        const { return grid_[idx];  }
    size_t NX() const { return grid_.NX(); }
    size_t NY() const { return grid_.NY(); }
    size_t NZ() const { return grid_.NZ(); }
    // -------------------------------------------
    void SetElement(int x,int y,int z,double v) { grid_.setGrid(x,y,z,v);     }
    /// Type definition of iterator over grid elements.
    typedef Grid<double>::iterator iterator;
    iterator begin() { return grid_.begin(); }
    iterator end()   { return grid_.end();   }
    /// Increment grid bin corresponding to point by given value.
    inline long int Increment(Vec3 const&, double);
    inline long int Increment(const double*, double);
    /// Increment grid bin by given value.
    inline long int Increment(int,int,int,double);
    /// \return grid value at specified bin.
    double GridVal(int x,int y,int z)        const { return grid_.element(x,y,z);   }
    /// \return grid index
    long int CalcIndex(int i, int j, int k) const { return grid_.CalcIndex(i,j,k); }
  private:
    Grid<double> grid_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// DataSet_GridDbl::Increment()
long int DataSet_GridDbl::Increment(Vec3 const& xyz, double f) {
  int i,j,k;
  if (CalcBins(xyz[0],xyz[1],xyz[2],i,j,k))
    return grid_.incrementBy(i,j,k,f);
  return -1L; 
}
// DataSet_GridDbl::Increment()
long int DataSet_GridDbl::Increment(const double* xyz, double f) {
  int i,j,k;
  if (CalcBins(xyz[0],xyz[1],xyz[2],i,j,k))
    return grid_.incrementBy(i,j,k,f);
  return -1L;
}
// DataSet_GridDbl::Increment()
long int DataSet_GridDbl::Increment(int i, int j, int k, double f) {
  return grid_.incrementBy(i,j,k,f);
}
#endif
