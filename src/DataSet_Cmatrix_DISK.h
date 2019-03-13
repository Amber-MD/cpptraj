#ifndef INC_DATASET_CMATRIX_DISK_H
#define INC_DATASET_CMATRIX_DISK_H
#include "DataSet_Cmatrix.h"
#include "NC_Cmatrix.h"
/// Used to hold pairwise distance matrix for clustering on disk.
class DataSet_Cmatrix_DISK : public DataSet_Cmatrix {
  public:
    DataSet_Cmatrix_DISK() : DataSet_Cmatrix(CMATRIX_DISK) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_Cmatrix_DISK(); }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return file_.MatrixSize(); }
    void Info()                          const { return;                    }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const {}
    int Allocate(SizeArray const&) { return 0; }
    /// \return Size in bytes of set
    size_t MemUsageInBytes() const { return 0; }
    // ----- Cmatrix functions -------------------
    /// \return an element indexed by sievedFrames.
    inline double GetFdist(int, int) const;
    /// Set element at column/row to given value
    void SetElement(int x, int y, double val)  { file_.WriteCmatrixElement(x, y, val); }
    /// \return Actual number of elements in matrix
    size_t Nelements()                   const { return file_.MatrixSize();            }
    /// \return size used by matrix in bytes
    size_t DataSize()                    const { return sievedFrames_.DataSize();      }
    /// \return Actual number of rows in the matrix
    size_t Nrows()                       const { return file_.MatrixRows();            }
    /// \return Element at given index.
    double GetElement(unsigned int idx)  const { return file_.GetCmatrixElement(idx);  }
    /// \return true if matrix needs setup
    bool NeedsSetup()                    const { return (file_.MatrixSize() < 1);      }
    /// \return true if matrix needs calculation
    bool NeedsCalc()                     const { return true;                          }
    /// No more distances will be added; flush to disk
    void Complete()                            { file_.Sync();                         }
  protected:
    int AllocateCmatrix(size_t);
    int SetCdist(ClusterDist*) { return 0; }
  private:
    NC_Cmatrix file_;
};
// ----- Inline functions ------------------------------------------------------
double DataSet_Cmatrix_DISK::GetFdist(int x, int y) const {
  return file_.GetCmatrixElement( sievedFrames_.FrameToIdx(x),
                                  sievedFrames_.FrameToIdx(y) );
}
#endif
