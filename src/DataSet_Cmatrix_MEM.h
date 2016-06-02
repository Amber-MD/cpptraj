#ifndef INC_DATASET_CMATRIX_MEM_H
#define INC_DATASET_CMATRIX_MEM_H
#include "DataSet_Cmatrix.h"
#include "Matrix.h"
/// Used to hold pairwise distance matrix for clustering in memory.
class DataSet_Cmatrix_MEM : public DataSet_Cmatrix {
  public:
    DataSet_Cmatrix_MEM() : DataSet_Cmatrix(CMATRIX) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_Cmatrix_MEM(); }
    /// Access internal matrix pointer to interface with file IO
    float*       Ptr()                         { return Mat_.Ptr();         }
    float const* Ptr()                   const { return Mat_.Ptr();         }
    /// Clear matrix
    //void Clear();
    // ----- DataSet functions -------------------
    size_t Size()                        const { return Mat_.size();        }
    void Info()                          const { return;                    }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Allocate(SizeArray const&);
    // ----- Cmatrix functions -------------------
    /// \return an element indexed by sievedFrames.
    inline double GetFdist(int, int) const;
    /// Set element at row/column to given value
    inline void SetElement(int, int, double);
    /// \return Actual number of elements in matrix
    size_t Nelements()        const { return Mat_.size();               }
    /// \return size used by matrix in bytes
    size_t DataSize() const;
    /// \return Actual number of rows in the matrix
    size_t Nrows() const { return Mat_.Nrows(); }
    // -------------------------------------------
    /// Allocate ignore array for given # of original rows and sieve.
    //int SetupIgnore(size_t, std::vector<char> const&, int);
    // Iterators to internal matrix.
    typedef Matrix<float>::iterator const_iterator;
    const_iterator begin() const { return Mat_.begin(); }
    const_iterator end()   const { return Mat_.end();   }
  protected:
    int AllocateCmatrix(size_t);
  private:
    Matrix<float> Mat_;
};
// Inline functions
double DataSet_Cmatrix_MEM::GetFdist(int col, int row) const {
  // row and col are based on original; convert to reduced
  // FIXME: This assumes GetElement will never be called for a sieved frame.
  return Mat_.element(sievedFrames_.FrameToIdx(col), 
                      sievedFrames_.FrameToIdx(row));
}

void DataSet_Cmatrix_MEM::SetElement(int col, int row, double val) {
  Mat_.setElement(col, row, val);
}
#endif
