#ifndef INC_DATASET_GRIDFLT_H
#define INC_DATASET_GRIDFLT_H
#include "DataSet_3D.h"
#include "Grid.h"
/// Single-precision three-dimensional grid.
class DataSet_GridFlt : public DataSet_3D {
  public:
    DataSet_GridFlt() : DataSet_3D(GRID_FLT, TextFormat(TextFormat::DOUBLE, 12, 4)) {}
    DataSet_GridFlt(DataSet_GridFlt const&);
    float& operator[](size_t idx)              { return grid_[idx];          }
    static DataSet* Alloc()       { return (DataSet*)new DataSet_GridFlt();  }
    Grid<float> const& InternalGrid()    const { return grid_; }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return grid_.size();        }
#   ifdef MPI
    // FIXME: Currently just sums up. Should this be a separate Sync function?
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info()                          const { return; }
    void WriteBuffer(CpptrajFile&,SizeArray const&) const;
    size_t MemUsageInBytes() const { return grid_.DataSize(); }
    // ----- DataSet_3D functions ----------------
    int Allocate3D(size_t x, size_t y, size_t z)          { return grid_.resize(x,y,z); }
    double GetElement(size_t x, size_t y, size_t z) const { return (double)grid_.element(x,y,z); }
    double operator[](size_t idx)                   const { return (double)grid_[idx];  }
    size_t NX() const { return grid_.NX(); }
    size_t NY() const { return grid_.NY(); }
    size_t NZ() const { return grid_.NZ(); }
    /// \return grid index
    long int CalcIndex(size_t i, size_t j, size_t k) const { return grid_.CalcIndex(i,j,k); }
    void ReverseIndex(long int n, size_t& i, size_t& j, size_t& k) const
      { return grid_.ReverseIndex(n,i,j,k); }
    void UpdateVoxel(long int i, double val) { grid_[i] += (float)val; }
    // -------------------------------------------
    void SetElement(size_t x,size_t y,size_t z,float v) { grid_.setGrid(x,y,z,v);     }
    /// Type definition of iterator over grid elements.
    typedef Grid<float>::iterator iterator;
    iterator begin() { return grid_.begin(); }
    iterator end()   { return grid_.end();   }
    /// Increment grid bin corresponding to point by given value.
    inline long int Increment(Vec3 const&, float);
    inline long int Increment(const double*, float);
    /// Increment grid bin by given value.
    inline long int Increment(size_t,size_t,size_t,float);
    /// \return grid value at specified bin.
    float GridVal(size_t x,size_t y,size_t z)        const { return grid_.element(x,y,z);   }
  private:
    Grid<float> grid_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// DataSet_GridFlt::Increment()
long int DataSet_GridFlt::Increment(Vec3 const& xyz, float f) {
  size_t i,j,k;
  if (Bin().Calc(xyz[0],xyz[1],xyz[2],i,j,k))
    return grid_.incrementBy(i,j,k,f);
  return -1L; 
}
// DataSet_GridFlt::Increment()
long int DataSet_GridFlt::Increment(const double* xyz, float f) {
  size_t i,j,k;
  if (Bin().Calc(xyz[0],xyz[1],xyz[2],i,j,k))
    return grid_.incrementBy(i,j,k,f);
  return -1L;
}
// DataSet_GridFlt::Increment()
long int DataSet_GridFlt::Increment(size_t i, size_t j, size_t k, float f) {
  return grid_.incrementBy(i,j,k,f);
}
#endif
