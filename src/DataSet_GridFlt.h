#ifndef INC_DATASET_GRIDFLT_H
#define INC_DATASET_GRIDFLT_H
#include "DataSet_3D.h"
#include "Grid.h"
/// Single-precision three-dimensional grid.
class DataSet_GridFlt : public DataSet_3D {
  public:
    DataSet_GridFlt() : DataSet_3D(GRID_FLT, 12, 4) {}
    float& operator[](size_t idx)              { return grid_[idx];         }
    static DataSet* Alloc()       { return (DataSet*)new DataSet_GridFlt(); }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return grid_.size();       }
    int Sync()                                 { return 1;                  }
    void Info()                          const { return;                    }
    // ----- DataSet_3D functions ----------------
    int Allocate3D(size_t x,size_t y,size_t z) { return grid_.resize(x,y,z); }
    void Write3D(CpptrajFile&,int,int,int) const;
    double GetElement(size_t x,size_t y,size_t z) const { return (double)grid_.element(x,y,z); }
    size_t NX() const { return grid_.NX(); }
    size_t NY() const { return grid_.NY(); }
    size_t NZ() const { return grid_.NZ(); }
    // -------------------------------------------
    void Increment(size_t x,size_t y,size_t z) { grid_.increment(x,y,z); }
  private:
    Grid<float> grid_;
};
#endif
