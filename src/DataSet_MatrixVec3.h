#ifndef INC_DATASET_MATRIXVEC3_H
#define INC_DATASET_MATRIXVEC3_H
#include "DataSet_2D.h"
#include "Matrix.h"
#include "Vec3.h"
/// Two-dimensional matrix of XYZ vectors.
class DataSet_MatrixVec3 : public DataSet_2D {
  public:
    DataSet_MatrixVec3() : DataSet_2D(MATRIX_VEC3, 12, 4) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_MatrixVec3();    }
    Vec3& operator[](size_t idx)               { return mat_[idx];          }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return mat_.size();        }
    int Sync()                                 { return 1;                  }
    void Info()                          const { return;                    }
    // ----- DataSet_2D functions ----------------
    int Allocate2D(size_t x,size_t y)          { return mat_.resize(x,y);   }
    int AllocateHalf(size_t x)                 { return mat_.resize(x,0L);  }
    int AllocateTriangle(size_t x)             { return mat_.resize(0L,x);  }
    void Write2D(CpptrajFile&, int, int) const;
    double GetElement(size_t x,size_t y) const { return mat_.element(x,y).Magnitude2(); }
    size_t Nrows()                       const { return mat_.Nrows();       }
    size_t Ncols()                       const { return mat_.Ncols();       }
    double* MatrixArray()                const;
    // -------------------------------------------
    int AddElement(const Vec3& v)              { return mat_.addElement(v); }
    void SetElement(size_t x,size_t y,const Vec3& v){ mat_.setElement(x,y,v);}
  private:
    Matrix<Vec3> mat_;
};
#endif
