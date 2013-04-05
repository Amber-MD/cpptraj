#ifndef INC_DATASET_MATRIXFLT_H
#define INC_DATASET_MATRIXFLT_H
#include "DataSet_2D.h"
#include "Matrix.h"
/// Single-precision two-dimensional matrix.
class DataSet_MatrixFlt : public DataSet_2D {
  public:
    DataSet_MatrixFlt() : DataSet_2D(MATRIX_FLT, 12, 4) {}
    float& operator[](size_t idx)              { return mat_[idx];          }
    static DataSet* Alloc() { return (DataSet*)new DataSet_MatrixFlt();     }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return mat_.size();        }
    int Sync()                                 { return 1;                  }
    void Info()                          const { return;                    }
    // ----- DataSet_2D functions ----------------
    int Allocate2D(size_t x,size_t y)          { return mat_.resize(x,y);   }
    int AllocateHalf(size_t x)                 { return mat_.resize(x,0L);  }
    int AllocateTriangle(size_t x)             { return mat_.resize(0L,x);  }
    void Write2D(CpptrajFile&, int, int) const;
    double GetElement(size_t x,size_t y) const { return (float)mat_.element(x,y);}
    size_t Nrows()                       const { return mat_.Nrows();       }
    size_t Ncols()                       const { return mat_.Ncols();       }
    double* MatrixArray()                const;
    DataSet_2D::MType Kind()             const { return (DataSet_2D::MType)mat_.Type(); }
    DataSet_2D::MatrixType Type()        const { return DataSet_2D::NO_OP;  }
    // -------------------------------------------
    int AddElement(float d)                    { return mat_.addElement(d); }
    void SetElement(size_t x,size_t y,float d) { mat_.setElement(x,y,d);    }
  private:
    Matrix<float> mat_;
};
#endif
