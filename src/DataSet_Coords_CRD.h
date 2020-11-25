#ifndef INC_DATASET_COORDS_CRD_H
#define INC_DATASET_COORDS_CRD_H
#include "DataSet_Coords.h"
#include "CompactFrameArray.h"
/// Hold a reduced-precision array of Frames
class DataSet_Coords_CRD : public DataSet_Coords {
  public:
    DataSet_Coords_CRD();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Coords_CRD(); }
    // ----- DataSet functions -------------------
    size_t Size()        const { return frames_.MaxFrames(); }
    SizeArray DimSizes() const { return SizeArray(1, Size()); }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info() const { return; }
    void Add(size_t, const void*) {}
    int Allocate(SizeArray const&);
    /// \return Size in bytes of set
    size_t MemUsageInBytes() const { return frames_.SizeInBytes(); }
    int MemAlloc(SizeArray const&);
    void CopyBlock(size_t, const DataSet*, size_t, size_t);
    // ----- DataSet_Coords functions ------------
    /// Add a frame.
    void AddFrame(Frame const&);
    /// Get a frame at position.
    void GetFrame(int, Frame&);
    /// Get a frame at position corresponding to mask.
    inline void GetFrame(int idx, Frame& fIn, AtomMask const& mIn) {
      //fIn.SetFromCRD( coords_[idx], mIn, numCrd_, numBoxCrd_, cInfo_.HasVel() );
    }
    /// Set CRD at position with frame.
    inline void SetCRD(int idx, Frame const& fIn) {
      //coords_[idx] = fIn.ConvertToCRD(numBoxCrd_, cInfo_.HasVel());
    }
    /// Set topology and coordinate information
    int CoordsSetup(Topology const&, CoordinateInfo const&);
    // -------------------------------------------
    /// \return estimated size in bytes for given # of frames.
    //size_t EstSizeInBytes(size_t nframes) const { return sizeInBytes(nframes, Top().Natom(), numBoxCrd_); }
  private:
    //static size_t sizeInBytes(size_t, size_t, size_t);

    CompactFrameArray frames_;
    int framesToReserve_; ///< Frames to reserve set by Allocate() for use in CoordsSetup()
};
#endif
