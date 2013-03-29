#ifndef INC_DATASET_MATRIXDBL_H
#define INC_DATASET_MATRIXDBL_H
#include "DataSet_2D.h"
#include "Matrix.h"
/// Double-precision two-dimensional matrix.
class DataSet_MatrixDbl : public DataSet_2D {
  public:
    DataSet_MatrixDbl() : DataSet_2D(MATRIX_DBL, 12, 4) {}
    double& operator[](size_t idx)             { return mat_[idx];          }
    static DataSet* Alloc() { return (DataSet*)new DataSet_MatrixDbl();     }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return mat_.size();        }
    int Sync()                                 { return 1;                  }
    void Info()                          const { return;                    }
    // ----- DataSet_2D functions ----------------
    int Allocate2D(size_t x,size_t y)          { return mat_.resize(x,y);   }
    int AllocateHalf(size_t x)                 { return mat_.resize(x,0L);  }
    int AllocateTriangle(size_t x)             { return mat_.resize(0L,x);  }
    void Write2D(CpptrajFile&, int, int) const;
    double GetElement(size_t x,size_t y) const { return mat_.element(x,y);  }
    size_t Nrows()                       const { return mat_.Nrows();       }
    size_t Ncols()                       const { return mat_.Ncols();       }
    double* MatrixArray()                const;
    // -------------------------------------------
    int AddElement(double d)                   { return mat_.addElement(d); }
    void SetElement(size_t x,size_t y,double d){ mat_.setElement(x,y,d);    }
  private:
    Matrix<double> mat_;
};
#endif
