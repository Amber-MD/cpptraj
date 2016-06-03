#ifndef INC_DATASET_CMATRIX_NOMEM_H
#define INC_DATASET_CMATRIX_NOMEM_H
#include "DataSet_Cmatrix.h"
#include "ClusterDist.h"
/// Used when no pairwise distance matrix is to be stored; all calcs done on the fly.
class DataSet_Cmatrix_NOMEM : public DataSet_Cmatrix {
  public:
    DataSet_Cmatrix_NOMEM() : DataSet_Cmatrix(CMATRIX_NOMEM), cdist_(0) {}
    ~DataSet_Cmatrix_NOMEM();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Cmatrix_NOMEM(); }
  // ----- DataSet functions -------------------
    size_t Size()                        const { return 0; }
    void Info()                          const { return;   }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    int Allocate(SizeArray const&) { return 0; }
    // ----- Cmatrix functions -------------------
    /// \return an element indexed by sievedFrames.
    inline double GetFdist(int f1, int f2) const { return cdist_->FrameDist(f1, f2); }
    /// Set element at column/row to given value
    void SetElement(int, int, double) {}
    /// \return Actual number of elements in matrix //TODO put in Cmatrix
    inline size_t Nelements() const;
    /// \return size used by matrix in bytes
    size_t DataSize() const { return sievedFrames_.DataSize(); }
    /// \return Actual number of rows in the matrix // TODO put in Cmatrix
    size_t Nrows() const { return sievedFrames_.ActualNframes(); }
    /// \return Element at given index.
    inline double GetElement(unsigned int) const;
    /// \return true if matrix needs setup
    bool NeedsSetup() const { return (cdist_ == 0); }
    // -------------------------------------------
    int SetupSieveAndCdist(size_t, size_t, int, ClusterDist*);
  protected:
    int AllocateCmatrix(size_t) { return 0; }
  private:
    ClusterDist* cdist_;
    ClusterSieve::SievedFrames idxToFrame_; // Needed for GetElement
};
// ----- Inline functions ------------------------------------------------------

size_t DataSet_Cmatrix_NOMEM::Nelements() const {
  return (sievedFrames_.ActualNframes()*(sievedFrames_.ActualNframes()-1)) / 2;
}

double DataSet_Cmatrix_NOMEM::GetElement(unsigned int idxIn) const {
  int iidx = (int)idxIn;
  return cdist_->FrameDist( idxToFrame_[ iidx / sievedFrames_.ActualNframes() ],
                            idxToFrame_[ iidx % sievedFrames_.ActualNframes() ] );
}
#endif
