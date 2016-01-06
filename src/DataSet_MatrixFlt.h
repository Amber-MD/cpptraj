#ifndef INC_DATASET_MATRIXFLT_H
#define INC_DATASET_MATRIXFLT_H
#include "DataSet_2D.h"
#include "Matrix.h"
/// Single-precision two-dimensional matrix.
class DataSet_MatrixFlt : public DataSet_2D {
  public:
    DataSet_MatrixFlt() : DataSet_2D(MATRIX_FLT, TextFormat(TextFormat::DOUBLE, 12, 4)) {}
    float& operator[](size_t idx)              { return mat_[idx];          }
    static DataSet* Alloc() { return (DataSet*)new DataSet_MatrixFlt();     }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return mat_.size();        }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    void Info()                          const { return;                    }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    // ----- DataSet_2D functions ----------------
    int Allocate2D(size_t x,size_t y)          { kind_=FULL; return mat_.resize(x,y); }
    int AllocateHalf(size_t x)                 { kind_=HALF; return mat_.resize(x,0); }
    int AllocateTriangle(size_t x)             { kind_=TRI;  return mat_.resize(0,x); }
    double GetElement(size_t x,size_t y) const { return (double)mat_.element(x,y);}
    double GetElement(size_t i)          const { return (double)mat_[i];    }
    size_t Nrows()                       const { return mat_.Nrows();       }
    size_t Ncols()                       const { return mat_.Ncols();       }
    double* MatrixArray()                const;
    MatrixKindType MatrixKind()          const { return kind_;              }
    // -------------------------------------------
    int AddElement(float d)                    { return mat_.addElement(d); }
    void SetElement(size_t x,size_t y,float d) { mat_.setElement(x,y,d);    }
  private:
    Matrix<float> mat_;
    MatrixKindType kind_;
};
#endif
