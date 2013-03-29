#ifndef INC_DATASET_MATRIXFLT_H
#define INC_DATASET_MATRIXFLT_H
#include "DataSet.h"
#include "DataSet_2D.h"
#include "Matrix.h"
/// Single-precision two-dimensional matrix.
class DataSet_MatrixFlt : public DataSet, DataSet_2D {
  public:
    DataSet_MatrixFlt() : DataSet(MATRIX2D, 12, 4, 2) {}
    float& operator[](size_t idx)              { return mat_[idx];          }
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
    // -------------------------------------------
    int AddElement(float d)                    { return mat_.addElement(d); }
    void SetElement(size_t x,size_t y,float d) { mat_.setElement(x,y,d);    }
    float* Ptr()                               { return mat_.Ptr();         }
  private:
    Matrix<float> mat_;
};
#endif
