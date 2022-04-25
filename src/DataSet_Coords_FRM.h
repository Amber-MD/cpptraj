#ifndef INC_DATASET_COORDS_FRM_H
#define INC_DATASET_COORDS_FRM_H
#include "DataSet_Coords.h"
/// Hold full-precision array of Frames in memory.
class DataSet_Coords_FRM : public DataSet_Coords {
  public:
    /// CONSTRUCTOR
    DataSet_Coords_FRM();
    /// Allocate DataSet
    static DataSet* Alloc() { return (DataSet*)new DataSet_Coords_FRM(); }
    // ----- DataSet functions -------------------
    /// \return Number of frames
    size_t Size()        const { return frames_.size(); }
    /// \return array containing sizes for each dimension (frames)
    SizeArray DimSizes() const { return SizeArray(1, Size()); }
#   ifdef MPI
    /// Synchronize all data to the master process
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    /// Print info to stdout
    void Info() const { return; }
    /// Add single element
    void Add(size_t, const void*);
    /// Reserve memory using given array of dimensions
    int Allocate(SizeArray const&);
    /// \return Size of set in bytes
    size_t MemUsageInBytes() const;
    /// Allocate array using given array of dimensions
    int MemAlloc(SizeArray const&);
    /// Copy a block of the DataSet
    void CopyBlock(size_t, const DataSet*, size_t, size_t);
    // ----- DataSet_Coords functions ------------
    /// Add a frame
    void AddFrame(Frame const&);
    /// Set frame at specified position
    void SetCRD(int, Frame const&);
    /// Get a frame at position.
    void GetFrame(int, Frame&);
    /// Get a frame at position corresponding to mask.
    void GetFrame(int, Frame&, AtomMask const&);
    /// Set topology and coordinate information
    int CoordsSetup(Topology const&, CoordinateInfo const&);
    // -------------------------------------------
  private:
    typedef std::vector<Frame> FrmArrayType;
    FrmArrayType frames_; ///< Array of Frames
};
#endif
