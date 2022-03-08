#ifndef INC_CLUSTER_CMATRIX_NC_H
#define INC_CLUSTER_CMATRIX_NC_H
#include <vector>
#include <string>
class FileName;
namespace Cpptraj {
namespace Cluster {
class Cframes;
/// NetCDF cluster matrix file.
class Cmatrix_NC {
  public:
    enum ModeType { READ=0, WRITE };
    Cmatrix_NC();
    ~Cmatrix_NC();
    /// \return Estimated disk usage in bytes for given # rows
    static size_t EstimatedDiskUsageInBytes(size_t);
    /// \return true if file is NetCDF cluster matrix file.
    static bool ID_Cmatrix(FileName const&);

#   ifdef BINTRAJ
    /// Open cluster matrix file for reading. Set sieve ID.
    int OpenCmatrixRead(FileName const&, int&);
    /// Get cluster matrix element (col, row)
    double GetCmatrixElement(unsigned int, unsigned int) const;
    /// Get cluster matrix element (raw index)
    double GetCmatrixElement(unsigned int) const;
    /// \return Array containing sieve status for each frame: 'T' sieved out, 'F' present
    std::vector<char> GetSieveStatus() const;
    /// Read cmatrix into given pointer
    int GetCmatrix(float*) const;
    /// Create cluster matrix file; # frames, # rows, sieve
    int CreateCmatrix(FileName const&, unsigned int, unsigned int, int, std::string const&);
    /// Sync to disk.
    void Sync() const;
    /// Reopen in shared write mode for random access
    int ReopenSharedWrite(FileName const&);
    /// Write non-sieved frames array.
    int WriteFramesArray(std::vector<int> const&) const; // TODO deprecate
    int WriteFramesArray(Cframes const&) const;
    /// Write cluster matrix element (col, row)
    int WriteCmatrixElement(unsigned int, unsigned int, double) const;
    /// Write cluster matrix using given pointer
    int WriteCmatrix(const float*) const;
    /// Close cluster matrix file.
    void CloseCmatrix();

    /// \return Matrix size
    unsigned int MatrixSize() const { return mSize_; }
    /// \return Matrix rows.
    unsigned int MatrixRows() const { return nRows_; }
    /// \return Current access mode
    ModeType Mode()           const { return mode_;  }
#   endif /* BINTRAJ */
  private:
#   ifdef BINTRAJ
    static inline bool IsCpptrajCmatrix(int);
    long int CalcIndex(unsigned int, unsigned int) const;

    int ncid_;                  ///< NetCDF file ID.
    int n_original_frames_DID_; ///< Number of original frames dimension.
    int n_rows_DID_;            ///< Number of rows (actual frames, N) in matrix dimension.
    int msize_DID_;             ///< Actual matrix size ( (N*(N-1))/2 ).
    int cmatrix_VID_;           ///< Cluster matrix variable ID ( matrix size ).
    int actualFrames_VID_;      ///< Non-sieved frames array ( N ).
    unsigned int nFrames_;      ///< Original number of frames.
    unsigned int nRows_;        ///< Number of rows (actual frames, N) in matrix dimension.
    unsigned int mSize_;        ///< Actual matrix size.
    ModeType mode_;             ///< Access mode.
#   endif /* BINTRAJ */
};

}
}
#endif
