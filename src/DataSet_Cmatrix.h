#ifndef INC_DATASET_CMATRIX_H
#define INC_DATASET_CMATRIX_H
#include "DataSet.h"
#include "ClusterSieve.h"
#include "ClusterDist.h"
/// Base class for pairwise distance matrices for clustering.
/** This matrix allows for sieving, i.e. it may hold less actual data than
  * the original size would warrant. For example, if there were originally
  * 10 rows of data the matrix would contain (10*9)/2 = 45 elements. However,
  * if every other frame is sieved (i.e. sieve 2) we really only need info
  * for (5*4)/2 = 10 elements; in this case we need to map original indices
  * to actual indices, which is what the ClusterSieve class is for.
  */
class DataSet_Cmatrix : public DataSet {
  public:
    DataSet_Cmatrix() {}
    DataSet_Cmatrix(DataType t) : DataSet(t, CLUSTERMATRIX,
                                          TextFormat(TextFormat::DOUBLE, 12, 4), 2) {}
    virtual ~DataSet_Cmatrix() {}
    // ----- DataSet functions -------------------
    // NOTE: Disabled for all DataSet_Cmatrix sets
    void Add(size_t, const void*) {}
    int Append(DataSet*) { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // ----- Cmatrix functions -------------------
    /// \return Distance between given frames (indexed by sievedFrames).
    virtual double GetFdist(int, int) const = 0;
    /// Set element at column/row to given value
    virtual void SetElement(int, int, double) = 0;
    /// \return Actual number of elements in matrix
    virtual size_t Nelements() const = 0;
    /// \return size used by matrix in bytes
    virtual size_t DataSize() const = 0;
    /// \return Actual number of rows in the matrix.
    virtual size_t Nrows() const = 0;
    /// \return Element at given index.
    virtual double GetElement(unsigned int) const = 0;
    /// \return true if matrix needs setup
    virtual bool NeedsSetup() const = 0;
    /// \return true if matrix needs calculation
    virtual bool NeedsCalc() const = 0;
    /// Indicate that no more distances will be added to matrix.
    virtual void Complete() = 0;
    // ----- Sieved frames functions -------------
    /// \return An array containing frame numbers that have not been sieved out.
    ClusterSieve::SievedFrames const& FramesToCluster() const { return sievedFrames_.Frames(); }
    /// \return Sieve value
    int SieveValue()                    const { return sievedFrames_.Sieve();  }
    /// \return Sieve type
    ClusterSieve::SieveType SieveType() const { return sievedFrames_.Type();   }
    /// \return Original number of frames before sieving.
    int OriginalNframes()               const { return sievedFrames_.MaxFrames(); }
    /// \return true if frame was sieved out.
    bool FrameWasSieved(int f)          const { return (sievedFrames_.FrameToIdx(f) == -1); }
    // -------------------------------------------
    /// Allocate matrix and sieve info for given size and sieve
    int SetupWithSieve(ClusterDist*, size_t, size_t, int);
    /// Allocate sieve info from given array.
    int SetSieveFromArray(std::vector<char> const&, int);
    /// Print matrix elements to STDOUT
    void PrintElements() const;
  protected:
    virtual int AllocateCmatrix(size_t) = 0;
    virtual int SetCdist(ClusterDist*) = 0;
    ClusterSieve sievedFrames_; ///< Hold info on frames actually being processed. TODO make private
};
#endif
